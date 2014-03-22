import re
import os
import subprocess
import logging

# Matches 1+ occurrences of a number, followed by a letter from {MIDNSHPX=}
CIGAR_RE = re.compile('[0-9]+[MIDNSHPX=]')

LOGGER = logging.getLogger(__name__)

# TODO:  handle X, =, P, N
def apply_cigar (cigar, seq, qual):
    """
    Parse SAM CIGAR and apply to the SAM nucleotide sequence.

    Input: cigar, sequence, and quality string from SAM.
    Output: shift (?), sequence with CIGAR incorporated + new quality string
    """

    newseq = ''
    newqual = ''
    tokens = CIGAR_RE.findall(cigar)
    if len(tokens) == 0:
        return None, None, None

    # Account for removing soft clipped bases on left
    shift = 0
    if tokens[0].endswith('S'):
        shift = int(tokens[0][:-1])

    left = 0
    for token in tokens:
        length = int(token[:-1])

        # Matching sequence: carry it over
        if token[-1] == 'M':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length

        # Deletion relative to reference: pad with gaps
        elif token[-1] == 'D':
            newseq += '-'*length
            newqual += ' '*length 		# Assign fake placeholder score (Q=-1)

        # Insertion relative to reference:
        elif token[-1] == 'I':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length
            continue

        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif token[-1] == 'S':
            left += length
            continue

        else:
            raise Exception("Unable to handle CIGAR token: {} - quitting".format(token))

    return shift, newseq, newqual


# TODO:  handle reverse complemented reads
# TODO:  handle reads where the insertsize is crappy
# TODO:  handle mate-pair vs paired-end
# TODO:  handle when reads align to multiple locations in genome
def merge_pairs (seq1, seq2, qual1, qual2, q_cutoff=10, minimum_q_delta=5):
    """
    Merge two sequences that overlap over some portion (paired-end
    reads).  Using the positional information in the SAM file, we will
    know where the sequences lie relative to one another.  In the case
    that the base in one read has no complement in the other read
    (in partial overlap region), take that base at face value.
    """

    mseq = ''

    # We swap the contents of the sequences so that seq2 is always the longest
    # so that when we iterate through the sequences, we don't run out of sequence before comparing the other
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1
        qual1, qual2 = qual2, qual1 # FIXME: quality strings must be concordant

    for i, c2 in enumerate(seq2):

        # FIXME: Track the q-score of each base at each position
        q2 = ord(qual2[i])-33

        if i < len(seq1):

            c1 = seq1[i]
            q1 = ord(qual1[i])-33

            # Gaps are given the lowest q-rank (q = -1) but are NOT N-censored
            if c1 == '-' and c2 == '-':
                mseq += '-'

            # NOT SURE IF THESE NEXT 3 CASES MAKE SENSE
            elif q1 is None and q2 is None:
                mseq += 'N'

            elif q2 is None and q1 is not None and q1 > q_cutoff:
                mseq += seq1[i]

            elif q1 is None and q2 is not None and q2 > q_cutoff:
                mseq += seq2[i]

            # Reads agree and one has sufficient confidence
            elif c1 == c2 and q1 > q_cutoff or q2 > q_cutoff:
                mseq += c1

            # Reads disagree but both have too similar a confidence
            elif abs(q2-q1) < minimum_q_delta:
                mseq += 'N'

            # Sequences disagree with differing confidence: take the high confidence
            elif q1 > q2 and q1 > q_cutoff:
                mseq += c1

            elif q2 > q1 and q2 > q_cutoff:
                mseq += c2

            # Not sure if I've missed any cases - for now, drop these cases....
            else:
                mseq += 'N'

        else:
            if q2 > q_cutoff:
                mseq += c2
            else:
                mseq += 'N'

    return mseq


def get_padded_seq_from_cigar(pos, cigar, seq, qual):
    """
    Returns the left-padded sequence from the cigar, with soft-clipped bases removed.
    Pads with '-' up to pos.

    :param int pos : pos field from SAM.  1-based position with respect to the reference.
    :param str cigar:  cigar field from SAM
    :param str seq : sequence field from SAM
    :param str qual : qual field from SAM
    :rtype str : left-padded sequence, with soft-clips removed

    """
    shift, formatted_seq, formatted_qual = apply_cigar(cigar, seq, qual)
    formatted_seq = '-' * (pos - 1) + formatted_seq
    formatted_qual = '-' * (pos - 1) + formatted_qual
    return [formatted_seq, formatted_qual]


def get_msa_fasta_from_sam(sam_filename, mapping_cutoff, read_qual_cutoff, max_prop_N, out_fasta_filename):
    """
    Parse SAM file contents for query-ref aligned sequences.
    Does pseudo multiple sequence alignment on all the query sequences and reference.
    TODO:  handle inserts.  Right now, all inserts are squelched so that there is multiple sequence alignment.
    For paired-end reads, merges the reads into a single sequence with gaps with respect to the reference.
    TODO:  handle mate pairs.
    Writes the MSA sequences to out_fasta_filename.

    :param sam_filename: full path to sam file
    :param mapping_cutoff:  Ignore alignments with mapping quality lower than the cutoff.
    :param read_qual_cutoff: When merging overlapping paired-end reads, ignore mate with read quality lower than the cutoff.
    :param max_prop_N:  Do not output sequences with proportion of N higher than the cutoff
    :param out_fasta_filename: full path of fasta file to write to.  Will completely overwite file.
    """

    LOGGER.debug("sam_filename=" + sam_filename + " out_fasta_filename=" + out_fasta_filename)
    with open(sam_filename, 'r') as sam_fh:
        with open(out_fasta_filename, 'w') as out_fasta_fh:
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

                # If read failed to map or has poor mapping quality, skip it
                ipos = int(pos)
                if refname == '*' or cigar == '*' or ipos == 0 or int(mapq) < mapping_cutoff:
                    continue

                padded_seq1, padded_qual1 = get_padded_seq_from_cigar(pos=ipos, cigar=cigar, seq=seq, qual=qual)

                padded_seq2 = ''
                padded_qual2 = ''
                if i < len(lines):
                    # Look ahead in the SAM for matching read
                    qname2, flag2, refname2, pos2, mapq2, cigar2, rnext2, pnext2, tlen2, seq2, qual2 = lines[i].rstrip().split('\t')[:11]

                    if qname2 == qname:  # TODO:  what if read maps multiple times?
                        i += 1

                        # If 2nd mate failed to map, then skip it
                        ipos2 = int(pos2)
                        if refname2 == '*' or cigar2 == '*' or ipos2 == 0 or int(mapq2) < mapping_cutoff:
                            continue

                        padded_seq2, padded_qual2 = get_padded_seq_from_cigar(pos=ipos2, cigar=cigar2, seq=seq2, qual=qual2)

                if padded_seq1 or padded_seq2:
                    # merge mates into one padded sequence
                    mseq = merge_pairs(padded_seq1, padded_seq2, padded_qual1, padded_qual2, read_qual_cutoff)

                    # Sequence must not have too many censored bases
                    if mseq.count('N') / float(len(mseq)) <= max_prop_N:
                        # Write multiple-sequence-aligned merged read to file using the name of the first mate
                        out_fasta_fh.write(">" + qname + "\n")
                        out_fasta_fh.write(mseq + "\n")



def samBitFlag(flag):
    """
    Interpret bitwise flag in SAM field as follows:

    Flag	Chr	Description
    =============================================================
    0x0001	p	the read is paired in sequencing
    0x0002	P	the read is mapped in a proper pair
    0x0004	u	the query sequence itself is unmapped
    0x0008	U	the mate is unmapped
    0x0010	r	strand of the query (1 for reverse)
    0x0020	R	strand of the mate
    0x0040	1	the read is the first read in a pair
    0x0080	2	the read is the second read in a pair
    0x0100	s	the alignment is not primary
    0x0200	f	the read fails platform/vendor quality checks
    0x0400	d	the read is either a PCR or an optical duplicate
    """
    labels = ['is_paired', 'is_mapped_in_proper_pair', 'is_unmapped', 'mate_is_unmapped',
              'is_reverse', 'mate_is_reverse', 'is_first', 'is_second', 'is_secondary_alignment',
              'is_failed', 'is_duplicate']

    binstr = bin(int(flag)).replace('0b', '')
    # flip the string
    binstr = binstr[::-1]
    # if binstr length is shorter than 11, pad the right with zeroes
    for i in range(len(binstr), 11):
        binstr += '0'

    bitflags = list(binstr)
    res = {}
    for i, bit in enumerate(bitflags):
        res.update({labels[i]: bool(int(bit))})

    return res


def create_depth_file_from_bam (bam_filename):
    """
    Gets the coverage from samtools depth.
    Creates a samtools per-base depth file with the same name as bam_filename but appended with ".depth".

    TODO: what to do with STDERR

    :param str bam_filename:  full path to sorted and indexed bam file
    :raise subprocess.CalledProcessError
    :rtype str : returns full filepath to depth file
    """

    import subprocess

    # Get per-base depth
    depth_filename = bam_filename + ".depth"
    subprocess.check_output(['samtools', 'depth', bam_filename], stdout=depth_fh, shell=False)
    return depth_filename


def sam_to_sort_bam(sam_filename, ref_filename):
    """
    Creates index of reference.  Converts sam to bam file sorted by coordinates.  Creates index of sorted bam.
    Pipes STDERR to STDOUT.
    Uses default samtools from PATH environment variable.

    :rtype str: full filepath to bam  file
    :param samtools_exe:  full filepath to samtools executable
    :param sam_filename: full filepath to sam alignment file
    """

    sam_filename_prefix, file_suffix = os.path.splitext(sam_filename)
    bam_filename = sam_filename_prefix + ".bam"

    # index reference
    index_ref_filename = ref_filename + ".fai"
    subprocess.check_call(['samtools', 'faidx', ref_filename], stderr=subprocess.STDOUT, shell=False)

    # convert sam to bam
    subprocess.check_call(['samtools', 'view', '-Sb', index_ref_filename, '-o', bam_filename, sam_filename],
                          stderr=subprocess.STDOUT, shell=False)

    # Sort the bam file by leftmost position on the reference assembly.  Required for samtools depth.
    sorted_bam_filename = bam_filename + ".sort"
    subprocess.check_call(['samtools', 'sort', '-f', bam_filename, sorted_bam_filename], stderr=subprocess.STDOUT, shell=False)

    # index sorted bam file
    subprocess.check_call(['samtools', 'index', sorted_bam_filename, ], stderr=subprocess.STDOUT, shell=False)




def get_ave_coverage_from_bam (bam_filename, ref, pos_start, pos_end):
    """
    If the depth file does not exist, then create one.
    Reads the depth file and returns average coverage in the specified region.

    :param str bam_filename: full path to bam file of alignments against the ref
    :param str ref: name of reference sequence to get coverage for
    :param int pos_start: region starting position, 1-based
    :param int pos_end: region ending position, 1-based
    :rtype float : average per-base coverage for the specified reference and region.
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
            if line_ref == ref and (line_pos >= pos_start or line_pos <= pos_end):
                total_coverage += line_depth

    ave_coverage = total_coverage / (pos_end - pos_start + 1)
    return ave_coverage
