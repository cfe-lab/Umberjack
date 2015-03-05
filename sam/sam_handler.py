import re
import os
import logging
from sam_constants import SamFlag as SamFlag
from sam_constants import  SamHeader as SamHeader
import sam_record
import paired_records
import Utility
import config.settings as settings      # sets the logging configs from logging.conf


NEWICK_NAME_RE = re.compile('[:;\-\(\)\[\]]')


LOGGER = logging.getLogger(__name__)


# TODO:  max_prop_N doesn't make sense if we include left and right pad gaps in the sequence length
# TODO:  how do we know that the sequence is on codon start????
def __write_seq(fh_out, name, seq, max_prop_N, breadth_thresh):
    """
    Helper function to write out sequence to fasta file handle if has sufficient bases.
    :param FileIO fh_out:  python file handle
    :param str name: Sequence Name
    :param str seq:  Sequence
    :param float max_prop_N: maximum fraction allowed N
    :param breadth_thresh:
    :return  True if sequence written out
    """
    # # Set stop codons to NNN so that HyPhy doesn't auto-remove a site without telling us
    # if is_mask_stop_codon:
    #     for nuc_pos in range(0, len(seq), Utility.NUC_PER_CODON):
    #         codon = seq[nuc_pos:nuc_pos+Utility.NUC_PER_CODON]
    #         if Utility.CODON2AA.get(codon, "") == Utility.STOP_AA:
    #             seq = seq[0:nuc_pos] + "NNN" + seq[nuc_pos+Utility.NUC_PER_CODON:]

    if seq.count('N') / float(len(seq)) <= max_prop_N and (seq.count("N") + seq.count("-"))/float(len(seq)) <= (1.0-breadth_thresh):
        # Newick tree formats don't like special characters.  Convert them to underscores.
        newick_nice_qname = re.sub(pattern=NEWICK_NAME_RE, repl='_', string=name)
        fh_out.write(">" + newick_nice_qname + "\n")
        fh_out.write(seq + "\n")
        return True
    return False


# TODO:  make stop codon removal optional
def create_msa_slice_from_sam(sam_filename, ref, out_fasta_filename, mapping_cutoff, read_qual_cutoff, max_prop_N,
                              breadth_thresh, start_pos=None, end_pos=None, is_insert=False, is_mask_stop_codon=False,
                              ref_len=None):
    """
    Parse SAM file contents for sequences aligned to a reference.
    Extracts the portion of the read that fits into the desired slice of the genome.
    For paired-end reads, merges the mates into a single sequence with gaps with respect to the reference.
    Creates a multiple sequence alignment (MSA) for the desired slice.
    Left and right pads the reads according to the positions within the slice.
    Writes the MSA sequences to out_fasta_filename.
    Converts query names so that they are compatible with Newick format in phylogenetic reconstruction by
        converting colons, semicolons, parentheses to underscores.



    :param is_mask_stop_codon:
    NB:  Sam file must be query sorted.
    NB:  Only takes the primary alignment.

    :param str sam_filename: full path to sam file.  Must be queryname sorted.
    :param str ref: name of reference contig to form MSA alignments to.
                    If None, then splits out all alignments to any reference that fit within the desired slice positions.
                    Setting the ref to None is only useful when the reads are aligned to a set of multiple sequence aligned
                    reference contigs, and you don't care which reference the read hits, just that it fits in the slice.
    :param str out_fasta_filename: full path to output multiple sequence aligned fasta file for the sequences in the slice.
    :param int mapping_cutoff:  Ignore alignments with mapping quality lower than the cutoff.
    :param int read_qual_cutoff: Convert bases with quality lower than this cutoff to N unless both mates agree.
    :param float max_prop_N:  Do not output reads with proportion of N higher than the cutoff.  Only counts the bases within the slice.
                                This only makes a difference if the slice is expected to be much wider than the (merged) read length.
    :param float breadth_thresh:  Fraction of the slice that the read must cover with actual bases A, C, G, T.
                                Reads below this threshold are excluded from output.
    :param int start_pos: 1-based start nucleotide start position of slice.  If None, then uses beginning of ref.
    :param int end_pos: 1-based end nucleotide start position of slice.  If None, then uses end of ref.
    :param bool is_insert: whether to exclude insertions to the reference.
                If include insertions, then the insertions will be multiple sequence aligned further by MAFFT.
    :param bool is_mask_stop_codon: whether to mask stop codons with "NNN".
                Most useful when you want to do codon analysis aftwards, as many codon models do not allow stop codons.
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
        for line in sam_fh:
            if line.startswith(SamHeader.TAG_HEADER_PREFIX):  # skip the headers
                continue

            lines_arr = line.rstrip().split('\t')
            if len(lines_arr) < 11:  # in case there are no alignments on this line
                continue

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
            #if (SamFlag.IS_UNMAPPED & flag or SamFlag.IS_SECONDARY_ALIGNMENT & flag or SamFlag.IS_CHIMERIC_ALIGNMENT & flag or
            #        rname == '*' or cigar == '*' or pos == 0 or mapq < mapping_cutoff or (ref and rname != ref)):
            if SamFlag.IS_UNMAPPED & flag or SamFlag.IS_SECONDARY_ALIGNMENT & flag or SamFlag.IS_CHIMERIC_ALIGNMENT & flag:
                continue

            is_mate_paired = SamFlag.IS_PAIRED & flag and not SamFlag.IS_MATE_UNMAPPED & flag

            mate = sam_record.SamRecord(ref_len=ref_len)
            mate.fill_record_vals(qname, flag, rname, seq, cigar, mapq, qual, pos, rnext, pnext)

            is_written = False

            if not prev_mate:  # We are not expecting this mate to be a pair with prev_mate.
                # If this record isn't paired and it passes our thresholds, then just write it out now
                if not is_mate_paired and mate.mapq >= mapping_cutoff and (not ref or mate.rname == ref):
                    # Do low quality masking on the this record and write it out.
                    mseq, mqual, stats = mate.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=True,
                                                           do_mask_low_qual=True, q_cutoff=read_qual_cutoff,
                                                           slice_start_wrt_ref_1based=start_pos,
                                                           slice_end_wrt_ref_1based=end_pos,
                                                           do_insert_wrt_ref=is_insert,
                                                           do_mask_stop_codon=is_mask_stop_codon)
                    is_written = __write_seq(out_fasta_fh, mate.qname, mseq, max_prop_N, breadth_thresh)
                    prev_mate = None
                # If this mate is paired, wait until we find the next mate in the pair or
                # know that the next mate doesn't exist before writing this mate out
                elif is_mate_paired:
                    prev_mate = mate
            elif prev_mate:  # We are expecting this mate to be a pair with prev_mate.
                if not SamFlag.IS_PAIRED & prev_mate.flag and not SamFlag.IS_MATE_UNMAPPED & prev_mate.flag:
                    raise ValueError("Invalid logic.  Shouldn't get here. - " +
                                     "Prevmate.qname=" + prev_mate.qname + " flag=" + str(prev_mate.flag) +
                                     " rnext=" + prev_mate.rnext +
                                     " mapq=" + str(prev_mate.mapq) + "\n" +
                                     "LINE=" + line)
                elif prev_mate.qname == mate.qname and not is_mate_paired:
                    raise ValueError("Previous mate and this mate have the same qname " + mate.qname +
                                     " but this mate isn't paired\n" +
                                     "LINE:" + line)
                elif prev_mate.qname == qname and is_mate_paired: # prev_mate and this mate are part of the same pair
                    # Only write out the mates with good map quality and hit our desired ref
                    if (prev_mate.mapq >= mapping_cutoff and mate.mapq >= mapping_cutoff and
                            (not ref or prev_mate.rname == ref or prev_mate.rname == "=") and
                            (not ref or mate.rname == ref or mate.rname == "=")):
                        pair = paired_records.PairedRecord(prev_mate, mate)
                        # TODO:  get a merge_sam_reads function that doesn't do all the stats
                        mseq, stats = pair.merge_sam_reads(q_cutoff=read_qual_cutoff,
                                                           pad_space_btn_mates="N",
                                                           do_insert_wrt_ref=is_insert,
                                                           do_pad_wrt_ref=False,
                                                           do_pad_wrt_slice=True,
                                                           do_mask_stop_codon=is_mask_stop_codon,
                                                           slice_start_wrt_ref_1based=start_pos,
                                                           slice_end_wrt_ref_1based=end_pos)
                        is_written = __write_seq(out_fasta_fh, pair.get_name(), mseq, max_prop_N, breadth_thresh)

                    elif prev_mate.mapq >= mapping_cutoff  and (not ref or prev_mate.rname == ref or prev_mate.rname == "="):
                        mseq, mqual, stats = prev_mate.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=True,
                                                                    do_mask_low_qual=True, q_cutoff=read_qual_cutoff,
                                                                    slice_start_wrt_ref_1based=start_pos,
                                                                    slice_end_wrt_ref_1based=end_pos,
                                                                    do_insert_wrt_ref=is_insert,
                                                                    do_mask_stop_codon=is_mask_stop_codon)
                        is_written = __write_seq(out_fasta_fh, prev_mate.qname, mseq, max_prop_N, breadth_thresh)

                    elif mate.mapq >= mapping_cutoff  and (not ref or mate.rname == ref or mate.rname == "="):
                        mseq, mqual, stats = mate.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=True,
                                                               do_mask_low_qual=True, q_cutoff=read_qual_cutoff,
                                                               slice_start_wrt_ref_1based=start_pos,
                                                               slice_end_wrt_ref_1based=end_pos,
                                                               do_insert_wrt_ref=is_insert,
                                                               do_mask_stop_codon=is_mask_stop_codon)
                        is_written = __write_seq(out_fasta_fh, mate.qname, mseq, max_prop_N, breadth_thresh)


                    # Clear the paired mate expectations for next set of sam records
                    prev_mate = None

                elif prev_mate.qname != mate.qname:  # This sam record does not pair with the previous sam record.
                    LOGGER.warn("Sam record inconsistent.  Expected pair for " + prev_mate.qname + " but got " + qname)

                    # Write out prev_mate
                    if prev_mate.mapq >= mapping_cutoff  and (not ref or prev_mate.rname == ref or prev_mate.rname == "="):
                        mseq, mqual, stats = prev_mate.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=True,
                                                                    do_mask_low_qual=True, q_cutoff=read_qual_cutoff,
                                                                    slice_start_wrt_ref_1based=start_pos,
                                                                    slice_end_wrt_ref_1based=end_pos,
                                                                    do_insert_wrt_ref=is_insert,
                                                                    do_mask_stop_codon=is_mask_stop_codon)
                        is_written = __write_seq(out_fasta_fh, prev_mate.qname, mseq, max_prop_N, breadth_thresh)

                    if is_mate_paired:  # May this mate will pair with the next record
                        prev_mate = mate
                else:
                    raise ValueError("Unpossible!")


            total_written += 1 if is_written else 0



        # In case the last record expects a mate that is not in the sam file
        # Do low quality masking on the previous record and write it out.
        if prev_mate:
            mseq, mqual = prev_mate.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=True, do_mask_low_qual=True,
                                                 q_cutoff=read_qual_cutoff, slice_start_wrt_ref_1based=start_pos,
                                                 slice_end_wrt_ref_1based=end_pos, do_insert_wrt_ref=is_insert,
                                                 do_mask_stop_codon=is_mask_stop_codon)
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
    """
    Searches sam header for reference length
    :param str sam_filename:  path to sam file
    :param str ref: name of reference
    :return: length of reference or None if not found in the header
    :rtype: int
    """
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
