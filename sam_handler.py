import re
import os
import subprocess
import logging


# Matches 1+ occurrences of a number, followed by a letter from {MIDNSHPX=}
CIGAR_RE = re.compile('[0-9]+[MIDNSHPX=]')
SEQ_PAD_CHAR = '-'
QUAL_PAD_CHAR = '!'     # This is the ASCII character for the lowest PHRED quality score in Sanger qualities
NEWICK_NAME_RE = re.compile('[:;\-\(\)\[\]]')
NUCL_RE = re.compile('[^nN\-]')
LOGGER = logging.getLogger(__name__)

PHRED_SANGER_OFFSET = 33

class SamFlag:
    IS_PAIRED =                0x001
    IS_MAPPED_IN_PROPER_PAIR = 0x002
    IS_UNMAPPED =              0x004
    IS_MATE_UNMAPPED =         0x008
    IS_REVERSE =               0x010
    IS_MATE_REVERSE =          0x020
    IS_FIRST =                 0x040
    IS_SECOND =                0x080
    IS_SECONDARY_ALIGNMENT =   0x100
    IS_FAILED =                0x200
    IS_DUPLICATE =             0x400
    IS_CHIMERIC_ALIGNMENT =    0x800


# TODO:  handle X, =, P, N
# TODO:  handle combinations of inserts/deletions
def apply_cigar (cigar, seq, qual):
    """
    Parse SAM CIGAR and apply to the SAM nucleotide sequence.
    Remove soft-clipped sequences.  They may be valid polymorphisms but there is no alignment information for clipped
    sequences so we have no way to know how they line up with other sequences.
    Bases with low quality are not removed here - that can be done with merge_pairs()
    Left-pads and right-pads sequence so that it lines up with the reference.

    :return:  tuple [left and right padded sequence, left and right padded quality]
    :rtype : tuple [str, str]
    :param str cigar: SAM cigar field
    :param str seq: SAM sequence field
    :param str qual: SAM quality field
    """

    newseq = ''
    newqual = ''
    tokens = CIGAR_RE.findall(cigar)
    if len(tokens) == 0:
        return None, None

    left = 0

    for token in tokens:
        length = int(token[:-1])

        if token[-1] == 'S':
            left += length

        # Matching sequence: carry it over
        elif token[-1] == 'M':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length

        # Deletion relative to reference: pad with gaps
        elif token[-1] == 'D':
            newseq += '-'*length
            newqual += '!'*length 		# Assign fake placeholder score (Q=-1)

        # Insertion relative to reference:
        elif token[-1] == 'I':
            # newseq += seq[left:(left+length)]
            # newqual += qual[left:(left+length)]
            left += length
            continue

        else:
            raise Exception("Unable to handle CIGAR token: {} - quitting".format(token))

    return newseq, newqual


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

    for i in range(0, len(seq2)):

        q2 = ord(qual2[i])-PHRED_SANGER_OFFSET
        q1 = ord(qual1[i])-PHRED_SANGER_OFFSET

        if q1 is None and q2 is None:
            mseq += 'N'

        elif q2 is None and q1 is not None and q1 > q_cutoff:
            mseq += seq1[i]

        elif q1 is None and q2 is not None and q2 > q_cutoff:
            mseq += seq2[i]

        # Reads agree  - quality doesn't matter
        elif seq1[i] == seq2[i]:  # Gaps are NOT N-censored
            mseq += seq1[i]

        # Sequences disagree with differing confidence: take the high confidence
        elif q1 > q2 and q1 > q_cutoff:
            mseq += seq1[i]

        elif q2 > q1 and q2 > q_cutoff:
            mseq += seq2[i]

        else:
            mseq += 'N'

    return mseq


# TODO:  handle inserts
# TODO:  hack - We hard cut sequences if they extend past the reference boundaries.  Don't do this.
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
    :return:  tuple [left-padded sequence with soft-clips removed, left-padded quality sequence with soft-clips removed]
    :rtype : tuple [str, str]

    """

    formatted_seq, formatted_qual = apply_cigar(cigar, seq, qual)
    left_pad_len = pos  - 1
    right_pad_len = ref_len - left_pad_len - len(formatted_seq)

    # TODO:  This hack is in place so that we don't have to worry about MSA alignments
    if right_pad_len < 0:
        formatted_seq = formatted_seq[:right_pad_len]
        formatted_qual = formatted_qual[:right_pad_len]
        right_pad_len = 0

    padded_seq = (SEQ_PAD_CHAR * left_pad_len) + formatted_seq + (SEQ_PAD_CHAR * right_pad_len)
    padded_qual = (QUAL_PAD_CHAR * left_pad_len) + formatted_qual + (QUAL_PAD_CHAR*right_pad_len)

    if len(padded_seq) != ref_len:
        raise Exception("len(padded_seq)=" + str(len(padded_seq)) + " ref2len[rname]=" + str(ref_len))

    return [padded_seq, padded_qual]


# TODO:  handle inserts.  Right now, all inserts are squelched so that there is multiple sequence alignment.
def create_msa_fasta_from_sam(sam_filename, ref, ref_len, out_fasta_filename, mapping_cutoff, read_qual_cutoff,
                              max_prop_N):
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
    :param str ref: name of reference contig to form MSA alignments to
    :param int ref_len: length of reference contig in nucleotides
    :param str out_fasta_filename: full path of fasta file to write to.  Will completely overwite file.
    :param int mapping_cutoff:  Ignore alignments with mapping quality lower than the cutoff.
    :param int read_qual_cutoff: Convert bases with quality lower than this cutoff to N unless both mates agree.
    :param float max_prop_N:  Do not output merged sequences with proportion of N higher than the cutoff
    """

    LOGGER.debug("sam_filename=" + sam_filename + " out_fasta_filename=" + out_fasta_filename)
    with open(sam_filename, 'r') as sam_fh, open(out_fasta_filename, 'w') as out_fasta_fh:
        lines = sam_fh.readlines()

        if not lines:
            LOGGER.warn("Empty intput SAM file " + sam_filename)

        # Skip top SAM header lines
        for start, line in enumerate(lines):
            if not line.startswith('@'):
                break

        i = start  # Keep track of the current line so that we can come back after traversing forward for mates
        while i < len(lines):
            qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[i].rstrip().split('\t')[:11]
            i += 1

            if not refname == ref:
                continue

            # If read failed to map or has poor mapping quality, skip it
            # Be careful!  The refname is only set to '*' and pos is only set to '0'
            #   if all reads in the mate pair are unmapped.
            # If the mate is mapped, but this read is not, the refname and pos will be set to the mate's reference hit.
            # From SAM specs:
            # "For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read
            #   should have RNAME and POS identical to its mate"
            if (SamFlag.IS_UNMAPPED & int(flag) or SamFlag.IS_SECONDARY_ALIGNMENT & int(flag) or
                    refname == '*' or cigar == '*' or int(pos) == 0 or int(mapq) < mapping_cutoff):
                continue


            padded_seq1, padded_qual1 = get_padded_seq_from_cigar(pos=int(pos), cigar=cigar, seq=seq, qual=qual,
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
                    if (refname2 == ref and qname2 == qname and
                            not(SamFlag.IS_UNMAPPED & int(flag2) or SamFlag.IS_SECONDARY_ALIGNMENT & int(flag2) or
                                refname2 == '*' or cigar2 == '*' or int(pos2) == 0 or int(mapq2) < mapping_cutoff)):
                        padded_seq2, padded_qual2 = get_padded_seq_from_cigar(pos=int(pos2), cigar=cigar2, seq=seq2,
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