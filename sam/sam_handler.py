import re
import os
import sys
import subprocess
import logging
import math
from sam_constants import SamFlag as SamFlag
from sam_constants import  SamHeader as SamHeader
import sam_record
import paired_records
import Utility



# Matches 1+ occurrences of a number, followed by a letter from {MIDNSHPX=}
CIGAR_RE = re.compile('[0-9]+[MIDNSHPX=]')
SEQ_PAD_CHAR = '-'
QUAL_PAD_CHAR = ' '     # This is the ASCII character right blow lowest PHRED quality score in Sanger qualities  (-1)
NEWICK_NAME_RE = re.compile('[:;\-\(\)\[\]]')
NUCL_RE = re.compile('[^nN\-]')

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)


PHRED_SANGER_OFFSET = 33






def apply_cigar (cigar, seq, qual, ref_start_pos):
    """
    Use CIGAR string (Compact Idiosyncratic Gapped Alignment Report) in SAM data
    to apply soft clips, insertions, and deletions to the read sequence.
    Any insertions relative to the sample consensus sequence are discarded to
    enforce a strict pairwise alignment, and returned separately in a
    dict object.
    :return:  tuple (left-clip length,
                left and right padded sequence,
                left and right padded quality,
                insertion dict  {1-based pos wrt ref right before insertions : inserted bases, inserted qual})
    :rtype : tuple (int, str, str, dict {int : (seq, seq))
    :param str cigar: SAM cigar field
    :param str seq: SAM sequence field
    :param str qual: SAM quality field
    :param int ref_start_pos:  1-based start position with respect to the reference of the unclipped portion of the sequence
    """
    newseq = ''
    newqual = ''
    insertions = {}
    tokens = CIGAR_RE.findall(cigar)
    if len(tokens) == 0:
        return None, None, None, None
    # Account for removing soft clipped bases on left
    shift = 0
    if tokens[0].endswith('S'):
        shift = int(tokens[0][:-1])
    left = 0  # position wrt sequence
    pos_wrt_ref = ref_start_pos
    for token in tokens:
        length = int(token[:-1])
        # Matching sequence: carry it over
        if token[-1] == 'M' or token[-1] == 'X' or token[-1] == '=':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length
            pos_wrt_ref += length
        # Deletion relative to reference: pad with gaps
        elif token[-1] == 'D' or token[-1] == 'P' or token[-1] == 'N':
            newseq += SEQ_PAD_CHAR*length
            newqual += QUAL_PAD_CHAR*length  # Assign fake placeholder score (Q=-1)
        # Insertion relative to reference: skip it (excise it)
        elif token[-1] == 'I':
            insertions.update({pos_wrt_ref-1: (seq[left:(left+length)], qual[left:(left+length)])})
            left += length
        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif token[-1] == 'S':
            left += length
        elif token[-1] == 'H':
            pass  # Hard clip does not leave sequence in SAM and pos already starts at unclipped pos
        else:
            raise ValueError("Unable to handle CIGAR token: {} - quitting".format(token))


    return shift, newseq, newqual, insertions




# TODO:  handle when reads align to multiple locations in genome
def merge_pairs(seq1, seq2, qual1, qual2, q_cutoff=10):
    """
    Merge two sequences that overlap over some portion (paired-end
    reads).  Using the positional information in the SAM file, we will
    know where the sequences lie relative to one another.  In the case
    that the base in one read has no complement in the other read
    (in partial overlap region), take that base at face value.

    :return:  merged paired-end read
    :rtype : str
    :param str seq1: first mate sequence in forward direction, left and right padded to line up with reference
    :param str seq2: second mate sequence in forward direction, left and right padded to line up with reference
    :param str qual1: first mate quality in forward direction, left and right padded to line up with reference.
                        We expect this to be sanger phred-33 based ASCII characters.
    :param str qual2: second mate quality in forward direction, left and right padded to line up with reference.
                        We expect this to be sanger phred-33 based ASCII characters.
    :param int q_cutoff: quality cutoff below which a base is converted to N if there is no consensus between the mates.
    """

    mseq = ''

    if len(seq1) != len(seq2):
        raise Exception("Expect left and right padded mate sequences such that they line up with the reference.  " +
                            "mate1=" + seq1 + ", mate2=" + seq2)

    # If overlapping mate bases agree, either mate quality can be lower as long as the
    #   Probability of single base error at the quality cutoff >= Probability of both mates having same erroneous base
    # Quality = -10 log P(error),  P(error) = 1e-(Quality/10)  by Phred Score definition.
    # Assume that the probability of a wrong base is independent across mates, positions.
    # Say the correct base = G, the permutations of wrong overlapping bases = {AA, AC, AT, CA, CC, CT, TA, TC, TT}.
    # ==> There are 3 permutations where the overlapping bases are the same; There are 9 possible permutations of wrong overlapping bases.
    #
    # P(single base error at q_cutoff) >= P(both mates same base & both mates wrong base)
    # P(single base error at q_cutoff) >= P(both mates wrong base) P(both mates same base | both mates wrong base)
    # 1e-(q_cutoff/10) >= P(mate1 1-base error) P(mate2 1-base error) [3 matching overlaps / 9 possible permutations of overlapping wrong bases)
    # 1e-(q_cutoff/10) >= 1e-q1/10 *  1e-q2/10  * 1/3
    # 3e-(q_cutoff/10) >= 1e-q1/10 *  1e-q2/10
    # log3 - q_cutoff/10 >= -q1/10 - q2/10
    # -10log3 + qcutoff <= q1 + q2
    # qcutoff_overlap = -10log3 + qcutoff <= q1 + q2
    q_cutoff_overlap = -10 * math.log10(3) + q_cutoff

    for i in range(0, len(seq2)):

        q2 = ord(qual2[i])-PHRED_SANGER_OFFSET
        q1 = ord(qual1[i])-PHRED_SANGER_OFFSET

        # Both mates have gap at this position
        if seq1[i] == seq2[i] == SEQ_PAD_CHAR:
            mseq += seq1[i]

        # Both mates have the same base at this position
        elif seq1[i] == seq2[i] and seq1[i] != SEQ_PAD_CHAR:
            if q_cutoff_overlap <= q1 + q2:
                mseq += seq1[i]
            else:
                mseq += "N"

        # Both mates have disagreeing bases at this position: take the high confidence
        # If one of mates has a gap, the gap quality is always lower than the threshold
        elif q1 > q2 and q1 >= q_cutoff:
                mseq += seq1[i]

        elif q2 > q1 and q2 >= q_cutoff:
                mseq += seq2[i]

        elif q1 == q2 or (q1 < q_cutoff and q2 < q_cutoff):  # TODO:  put in the correct pairs capabilities into here
            mseq += "N"

        else:
            raise ValueError("We should never get here.  Unanticipated use case for merging sam reads")

    return mseq




# TODO:  handle inserts
# TODO:  hack - We hard cut sequences if they extend past the reference boundaries.  Don't do this.
# TODO:  ==> WTF?  sequences should never extend past ref!!!  when did this ever happen?
def get_padded_seq_from_cigar(pos, cigar, seq, qual, ref_len):
    """
    Returns the padded sequence from the cigar, with soft-clipped bases removed.
    Left Pads with '-' up to pos.
    Right pads with '-' until the end of the reference.

    :param int pos : pos field from SAM.  1-based position with respect to the reference.
    :param str cigar:  cigar field from SAM
    :param str seq : sequence field from SAM
    :param str qual : qual field from SAM
    :param int ref_len : length of reference contig in nucleotides
    :return:  tuple [left-padded sequence with soft-clips removed,
                left-padded quality sequence with soft-clips removed,
                insertion dict {1-based position with respect to reference right after inserts: (inserted seq, inserted qual)},
                left pad length,
                right pad length
    ]
    :rtype : tuple [str, str, dict {int:(str, str), int, int ]
    """

    left_clip_len, formatted_seq, formatted_qual, ref_pos_to_insert_seq_qual = apply_cigar(cigar, seq, qual, pos)
    left_pad_len = pos  - 1
    right_pad_len = ref_len - left_pad_len - len(formatted_seq)

    # TODO:  This hack is in place so that we don't have to worry about MSA alignments
    # TODO:  wtf?  we're limiting the seuqence to - right_pad_len bases before the end of the sequence.  Methinks this is a bug, but it never gets here.
    # if right_pad_len < 0:
    #     formatted_seq = formatted_seq[:right_pad_len]
    #     formatted_qual = formatted_qual[:right_pad_len]
    #     right_pad_len = 0

    padded_seq = (SEQ_PAD_CHAR * left_pad_len) + formatted_seq + (SEQ_PAD_CHAR * right_pad_len)
    padded_qual = (QUAL_PAD_CHAR * left_pad_len) + formatted_qual + (QUAL_PAD_CHAR*right_pad_len)

    if len(padded_seq) != ref_len:
        raise Exception("len(padded_seq)=" + str(len(padded_seq)) + " ref2len[rname]=" + str(ref_len))

    return padded_seq, padded_qual, ref_pos_to_insert_seq_qual, left_pad_len, right_pad_len


# TODO:  handle inserts.  Right now, all inserts are squelched so that there is multiple sequence alignment.
def create_msa_fasta_from_sam(sam_filename, ref, out_fasta_filename, mapping_cutoff, read_qual_cutoff,
                              max_prop_N, ref_len=None):
    """
    Parse SAM file contents for query-ref (pairwise) aligned sequences for a specific reference contig.
    For paired-end reads, merges the reads into a single sequence with gaps with respect to the reference.
    Creates a pseudo-multiple sequence alignment on all the query sequences and reference.
    Writes the MSA sequences to out_fasta_filename.
    Converts query names so that they are compatible with Newick format in phylogenetic reconstruction by
        converting colons, semicolons, parentheses to underscores.

    ASSUMES:  SAM is ordered by query name
    NB: Only takes the primary alignment.

    :param str sam_filename: full path to sam file
    :param str ref: name of reference contig to form MSA alignments to.  If None, then splits out all alignments to any reference.
    :param int ref_len: length of reference.  If None, then takes length from sam headers.
    :param str out_fasta_filename: full path of fasta file to write to.  Will completely overwite file.
    :param int mapping_cutoff:  Ignore alignments with mapping quality lower than the cutoff.
    :param int read_qual_cutoff: Convert bases with quality lower than this cutoff to N unless both mates agree.
    :param float max_prop_N:  Do not output merged sequences with proportion of N higher than the cutoff
    """

    if not ref_len:
        ref_len = get_ref_len(sam_filename, ref)

    LOGGER.debug("sam_filename=" + sam_filename + " out_fasta_filename=" + out_fasta_filename)
    with open(sam_filename, 'r') as sam_fh, open(out_fasta_filename, 'w') as out_fasta_fh:
        lines = sam_fh.readlines()

        if not lines:
            LOGGER.warn("Empty intput SAM file " + sam_filename)

        start = 0
        # Skip top SAM header lines
        for start, line in enumerate(lines):
            if not line.startswith('@'):
                break

        i = start  # Keep track of the current line so that we can come back after traversing forward for mates
        while i < len(lines):

            lines_arr = lines[i].rstrip().split('\t')
            if len(lines_arr) < 11:  # in case there are no alignments
                break

            qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines_arr[:11]
            i += 1

            if ref and not refname == ref:
                continue

            # If read failed to map or has poor mapping quality, skip it
            # Be careful!  The refname is only set to '*' and pos is only set to '0'
            #   if all reads in the mate pair are unmapped.
            # If the mate is mapped, but this read is not, the refname and pos will be set to the mate's reference hit.
            # From SAM specs:
            # "For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read
            #   should have RNAME and POS identical to its mate"
            if (SamFlag.IS_UNMAPPED & int(flag) or SamFlag.IS_SECONDARY_ALIGNMENT & int(flag) or
                    refname == '*' or cigar == '*' or int(pos) == 0 or mapq < mapping_cutoff):
                continue


            padded_seq1, padded_qual1, ref_pos_to_insert_seq_qual, left_pad_len, right_pad_len = get_padded_seq_from_cigar(pos=int(pos), cigar=cigar, seq=seq, qual=qual,
                                                                  ref_len=ref_len)

            padded_seq2 = ''
            padded_qual2 = ''
            if not SamFlag.IS_MATE_UNMAPPED & int(flag):
                #while i < len(lines) and padded_seq2 == '' and qname2 == qname:
                if i < len(lines):
                    # Look ahead in the SAM for matching read
                    qname2, flag2, refname2, pos2, mapq2, cigar2, rnext2, pnext2, tlen2, seq2, qual2 = lines[i].rstrip().split('\t')[:11]
                    i += 1

                    # if not refname2 == ref:
                    #     continue
                    #
                    # # If not the 2nd mate, failed to map, or is secondary alignment, then skip it
                    # if (not qname2 == qname or
                    #         SamFlag.IS_UNMAPPED & int(flag2) or SamFlag.IS_SECONDARY_ALIGNMENT & int(flag2) or
                    #         refname2 == '*' or cigar2 == '*' or int(pos2) == 0 or int(mapq2) < mapping_cutoff):
                    #     continue
                    #
                    # padded_seq2, padded_qual2 = get_padded_seq_from_cigar(pos=int(pos2), cigar=cigar2, seq=seq2,
                    #                                                       qual=qual2, flag=flag2, rname=refname2,
                    #                                                       ref_len=ref_len)
                    if ((not ref or refname2 == ref) and qname2 == qname and
                            not(SamFlag.IS_UNMAPPED & int(flag2) or SamFlag.IS_SECONDARY_ALIGNMENT & int(flag2) or
                                refname2 == '*' or cigar2 == '*' or int(pos2) == 0 or int(mapq2) < mapping_cutoff)):
                        padded_seq2, padded_qual2, ref_pos_to_insert_seq_qual, left_pad_len, right_pad_len = get_padded_seq_from_cigar(pos=int(pos2), cigar=cigar2, seq=seq2,
                                                                              qual=qual2, ref_len=ref_len)

            if padded_seq1 and padded_seq2:
                # merge mates into one padded sequence
                # We merge in case the mates overlap, the overlap gives us confidence on the bases
                mseq = merge_pairs(padded_seq1, padded_seq2, padded_qual1, padded_qual2, read_qual_cutoff)

                # Sequence must not have too many censored bases
                if mseq.count('N') / float(len(mseq)) <= max_prop_N:
                    # Write multiple-sequence-aligned merged read to file
                    # Newick tree formats don't like special characters.  Convert them to underscores.
                    newick_nice_qname = re.sub(pattern=NEWICK_NAME_RE, repl='_', string=qname)
                    out_fasta_fh.write(">" + newick_nice_qname + "\n")
                    out_fasta_fh.write(mseq + "\n")
            else:
                LOGGER.warn("Sequence has no mate " + qname)

# TODO:  make stop codon remove optional
def __write_seq(fh_out, name, seq, max_prop_N, breadth_thresh):
    """
    Helper function to write out sequence to fasta file handle if has sufficient bases.
    :return  True if sequence written out
    :rtype : boolean
    :param fh_out:
    :param name:
    :param seq:
    """
    # Set stop codons to NNN so that HyPhy doesn't auto-remove a site without telling us
    for nuc_pos in range(0, len(seq), Utility.NUC_PER_CODON):
        codon = seq[nuc_pos:nuc_pos+Utility.NUC_PER_CODON]
        if Utility.CODON2AA.get(codon, "") == Utility.STOP_AA:
            seq = seq[0:nuc_pos] + "NNN" + seq[nuc_pos+Utility.NUC_PER_CODON:]

    if seq.count('N') / float(len(seq)) <= max_prop_N and (seq.count("N") + seq.count("-"))/float(len(seq)) <= (1.0-breadth_thresh):
        # Newick tree formats don't like special characters.  Convert them to underscores.
        newick_nice_qname = re.sub(pattern=NEWICK_NAME_RE, repl='_', string=name)
        fh_out.write(">" + newick_nice_qname + "\n")
        fh_out.write(seq + "\n")
        return True
    return False


# TODO:  make stop codon removal optional
def create_msa_slice_from_sam(sam_filename, ref, out_fasta_filename, mapping_cutoff, read_qual_cutoff,
                              max_prop_N, breadth_thresh, start_pos=None, end_pos=None, is_insert=False, ref_len=None):
    """
    Parse SAM file contents for query-ref (pairwise) aligned sequences for a specific reference contig.
    For paired-end reads, merges the reads into a single sequence with gaps with respect to the reference.
    Creates a pseudo-multiple sequence alignment on all the query sequences and reference.
    Writes the MSA sequences to out_fasta_filename.
    Converts query names so that they are compatible with Newick format in phylogenetic reconstruction by
        converting colons, semicolons, parentheses to underscores.


    Sam file must be query sorted.  Secondary alignments are ignored.
    NB: Only takes the primary alignment.

    :param str sam_filename: full path to sam file.  Must be queryname sorted.
    :param str ref: name of reference contig to form MSA alignments to.  If None, then splits out all alignments to any reference.
    :param str out_fasta_filename: full path of fasta file to write to.  Will completely overwite file.
    :param int mapping_cutoff:  Ignore alignments with mapping quality lower than the cutoff.
    :param int read_qual_cutoff: Convert bases with quality lower than this cutoff to N unless both mates agree.
    :param float max_prop_N:  Do not output merged sequences with proportion of N higher than the cutoff
    :param int start_pos: 1-based start nucleotide start position of slice.  If None, then uses beginning of ref.
    :param int end_pos: 1-based end nucleotide start position of slice.  If None, then uses end of ref.
    :param int ref_len: length of reference.  If None, then takes length from sam headers.
    """


    LOGGER.debug("About to slice fasta " + out_fasta_filename + " from " + sam_filename)
    if os.path.exists(out_fasta_filename) and os.path.getsize(out_fasta_filename):
        LOGGER.warn("Found existing Sliced MSA-Fasta " + out_fasta_filename + ". Not regenerating.")
        total_seq = Utility.get_total_seq_from_fasta(out_fasta_filename)
        LOGGER.debug("Done slice fasta " + out_fasta_filename)
        return total_seq

    if not is_query_sort(sam_filename):
        raise ValueError("Sam file must be queryname sorted and header must specify sort order")

    total_written = 0
    if not ref_len:
        ref_len = get_reflen(sam_filename, ref)
    with open(sam_filename, 'r') as sam_fh, open(out_fasta_filename, 'w') as out_fasta_fh:
        prev_mate = None
        is_expect_mate = False
        for line in sam_fh:
            if line.startswith(SamHeader.TAG_HEADER_PREFIX):  # skip the headers
                continue

            lines_arr = line.rstrip().split('\t')
            if len(lines_arr) < 11:  # in case there are no alignments
                break

            qname, flag_str, rname, pos_str, mapq_str, cigar, rnext, pnext, tlen, seq, qual = lines_arr[:11]
            flag = int(flag_str)
            pos = int(pos_str)
            mapq = int(mapq_str)


            # If read failed to map or has poor mapping quality, skip it
            # Be careful!  The rname is only set to '*' and pos is only set to '0'
            #   if all reads in the mate pair are unmapped.
            # If the mate is mapped, but this read is not, the rname and pos will be set to the mate's reference hit.
            # From SAM specs:
            # "For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read
            #   should have RNAME and POS identical to its mate"
            if (SamFlag.IS_UNMAPPED & flag or SamFlag.IS_SECONDARY_ALIGNMENT & flag or SamFlag.IS_CHIMERIC_ALIGNMENT & flag or
                    rname == '*' or cigar == '*' or pos == 0 or mapq < mapping_cutoff or (ref and rname != ref)):
                continue

            mate = sam_record.SamRecord(ref_len=ref_len)
            mate.fill_record_vals(qname, flag, rname, seq, cigar, mapq, qual, pos, rnext, pnext)

            is_mate_paired = SamFlag.IS_PAIRED & flag and not SamFlag.IS_MATE_UNMAPPED & flag and (rnext == "=" or not ref or rnext == ref)

            # We are expecting this record or future records to be the next mate in pair
            if is_expect_mate:
                if (not prev_mate or not SamFlag.IS_PAIRED & prev_mate.flag or SamFlag.IS_MATE_UNMAPPED & prev_mate.flag or
                         not (prev_mate.rnext == "=" or not ref or prev_mate.rnext == ref)):
                    if prev_mate:
                        LOGGER.error("Invalid logic.  Prevmate.qname=" + prev_mate.qname + " flag=" + str(prev_mate.flag) + " rnext=" + prev_mate.rnext)
                    raise ValueError("Invalid logic.  Shouldn't get here. - " + line)

                # Check if last sam record and this sam record are for the same paired read
                if prev_mate.qname == qname and is_mate_paired:
                    pair = paired_records.PairedRecord(prev_mate, mate)
                    # TODO:  get a merge_sam_reads functioin that doesn't do all the stats
                    mseq, stats = pair.merge_sam_reads(q_cutoff=read_qual_cutoff,
                                                       pad_space_btn_mates="N",
                                                       do_insert_wrt_ref=is_insert,
                                                       do_pad_wrt_ref=False,
                                                       do_pad_wrt_slice=True,
                                                       slice_start_wrt_ref_1based=start_pos,
                                                       slice_end_wrt_ref_1based=end_pos)
                    is_written = __write_seq(out_fasta_fh, pair.get_name(), mseq, max_prop_N, breadth_thresh)
                    total_written += 1 if is_written else 0
                else:
                    # TODO:  dont spit out this warning if we are missing a mate because we rejected its map qual
                    LOGGER.warn("Sam record inconsistent.  Expected pair for " + prev_mate.qname + " but got " + qname)
                    # This sam record does not pair with the previous sam record.
                    # Do low quality masking on the previous record and write it out.
                    mseq, mqual, stats = prev_mate.get_seq_qual(do_pad_wrt_ref=False,
                                                         do_pad_wrt_slice=True,
                                                         do_mask_low_qual=True,
                                                         q_cutoff=read_qual_cutoff,
                                                         slice_start_wrt_ref_1based=start_pos,
                                                         slice_end_wrt_ref_1based=end_pos,
                                                         do_insert_wrt_ref=is_insert)
                    is_written = __write_seq(out_fasta_fh, prev_mate.qname, mseq, max_prop_N, breadth_thresh)
                    total_written += 1 if is_written else 0

            if not is_mate_paired:
                # Do low quality masking on the this record and write it out.
                mseq, mqual, stats = mate.get_seq_qual(do_pad_wrt_ref=False,
                                                do_pad_wrt_slice=True,
                                                do_mask_low_qual=True,
                                                q_cutoff=read_qual_cutoff,
                                                slice_start_wrt_ref_1based=start_pos,
                                                slice_end_wrt_ref_1based=end_pos,
                                                do_insert_wrt_ref=is_insert)
                is_written = __write_seq(out_fasta_fh, mate.qname, mseq, max_prop_N, breadth_thresh)
                total_written += 1 if is_written else 0
                is_expect_mate = False
                prev_mate = None
            else:
                if is_expect_mate and prev_mate.qname == qname and is_mate_paired:  # already wrote the pair
                    is_expect_mate = False
                    prev_mate = None
                else:
                    is_expect_mate = True
                    prev_mate = mate



        # In case we went through the full sam file but couldn't find the next paired mate to map well to the same ref
        # Do low quality masking on the previous record and write it out.
        if is_expect_mate and prev_mate:
            mseq, mqual = prev_mate.get_seq_qual(do_pad_wrt_ref=False,
                                                 do_pad_wrt_slice=True,
                                                 do_mask_low_qual=True,
                                                 q_cutoff=read_qual_cutoff,
                                                 slice_start_wrt_ref_1based=start_pos,
                                                 slice_end_wrt_ref_1based=end_pos,
                                                 do_insert_wrt_ref=is_insert)
            is_written = __write_seq(out_fasta_fh, prev_mate.qname, mseq, max_prop_N, breadth_thresh)
            total_written += 1 if is_written else 0

    LOGGER.debug("Done slice fasta " + out_fasta_filename)
    return total_written



def is_query_sort(sam_filename):
    """
    :param sam_filename:
    :return: whether the sam is sorted by queryname using the header.  False if header if absent.
    :rtype: boolean
    """
    with open(sam_filename, 'rU') as fh_in:
        header = fh_in.next()  # header should be the first line if present
        header = header.rstrip()
        if not header.startswith(SamHeader.TAG_HEADER_START):
            return False
        tags = header.split(SamHeader.TAG_SEP)
        for tag in tags[1:]:
            key, val = tag.split(SamHeader.TAG_KEY_VAL_SEP)
            if key == SamHeader.TAG_SORT_ORDER_KEY:
                if val == SamHeader.TAG_SORT_ORDER_VAL_QUERYNAME:
                    return True
                else:
                    return False
    return False

def get_reflen(sam_filename, ref):
    with open(sam_filename, 'rU') as fh_in:
        for line in fh_in:
            if not line.startswith(SamHeader.TAG_HEADER_PREFIX):
                return None
            if line.startswith(SamHeader.TAG_REFSEQ_START):
                line = line.rstrip()
                is_found_ref = False
                length = None
                for tag in line.split(SamHeader.TAG_SEP)[1:]:
                    key, val = tag.split(SamHeader.TAG_KEY_VAL_SEP)
                    if key == SamHeader.TAG_REFSEQ_NAME_KEY and val == ref:
                        is_found_ref = True
                    if key == SamHeader.TAG_REFSEQ_LEN_KEY:
                        length = int(val)
                if is_found_ref:
                    return length
    return None



def create_depth_file_from_bam(bam_filename):
    """
    Gets the coverage from samtools depth.
    Creates a samtools per-base depth file with the same name as bam_filename but appended with ".depth".

    TODO: what to do with STDERR

    :param str bam_filename:  full path to sorted and indexed bam file
    :return: returns full filepath to depth file
    :rtype : str
    :raise subprocess.CalledProcessError
    """

    # Get per-base depth
    depth_filename = bam_filename + ".depth"
    with open(depth_filename, 'w') as depth_fh:
        subprocess.check_call(['samtools', 'depth', bam_filename], stdout=depth_fh, shell=False)
    return depth_filename


def sam_to_sort_bam(sam_filename, ref_filename):
    """
    Creates index of reference.  This creates a <ref_filename>.fai file.
    Converts sam to bam file sorted by coordinates.
    This creates a <sam_filename prefix>.bam and <sam filename prefix>.bam.sort files.
    Creates index of sorted bam.  This creates a <bam_filename>.index file.
    Pipes STDERR to STDOUT.
    Uses default samtools from PATH environment variable.

    :return: full filepath to sorted, indexed bam  file
    :rtype : str
    :param str sam_filename: full filepath to sam alignment file
    :param str ref_filename:  full filepath to reference fasta
    """

    sam_filename_prefix, file_suffix = os.path.splitext(sam_filename)
    bam_filename = sam_filename_prefix + ".bam"

    # index reference
    index_ref_filename = ref_filename + ".fai"
    subprocess.check_call(['samtools', 'faidx', ref_filename], stderr=subprocess.STDOUT, shell=False)

    # convert sam to bam
    subprocess.check_call(['samtools', 'view', '-Sb', '-o', bam_filename, sam_filename],
                          stderr=subprocess.STDOUT, shell=False)

    # Sort the bam file by leftmost position on the reference assembly.  Required for samtools depth.
    sorted_bam_filename_prefix = sam_filename_prefix + ".sort"
    sorted_bam_filename = sorted_bam_filename_prefix + ".bam"
    subprocess.check_call(['samtools', 'sort', bam_filename, sorted_bam_filename_prefix], stderr=subprocess.STDOUT, shell=False)

    # index sorted bam file
    subprocess.check_call(['samtools', 'index', sorted_bam_filename, ], stderr=subprocess.STDOUT, shell=False)

    return sorted_bam_filename


def get_ave_coverage_from_bam (bam_filename, ref, pos_start, pos_end):
    """
    If the depth file does not exist, then create one.
    Reads the depth file and returns average coverage in the specified region.

    :param str bam_filename: full path to bam file of alignments against the ref
    :param str ref: name of reference sequence to get coverage for
    :param int pos_start: region starting position, 1-based
    :param int pos_end: region ending position, 1-based
    :return: average per-base coverage for the specified reference and region.
    :rtype : float
    """
    import os

    depth_filename = bam_filename + ".depth"
    if not os.path.isfile(depth_filename) or os.path.getsize(depth_filename) <= 0:
        create_depth_file_from_bam(bam_filename)

    total_coverage = 0
    with open(depth_filename, 'r') as depth_fh:
        for line in depth_fh:
            # Output of samtools depth:
            # <reference>   <position>  <depth>
            line_ref, line_pos, line_depth = line.rstrip().split('\t')
            if line_ref == ref and int(line_pos) >= pos_start and int(line_pos) <= pos_end:
                total_coverage += int(line_depth)

    ave_coverage = round(float(total_coverage) / (pos_end - pos_start + 1), 0)
    return ave_coverage


def get_ref_len(sam_filename, ref):
    with open(sam_filename, 'rU') as fh:
        for line in fh:
            if line.startswith("@SQ"):
                fields = line.split("\t")
                if line.find("SN:" + ref) >= 0:
                    line_reflen_match = re.search(r"\tLN:(\d+)", line)
                    return int(line_reflen_match.group(1))

    return None