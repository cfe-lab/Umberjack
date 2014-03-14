import re

# Matches 1+ occurrences of a number, followed by a letter from {MIDNSHPX=}
CIGAR_RE = re.compile('[0-9]+[MIDNSHPX=]')

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


def get_seq_from_sam(sam_filename, mapping_cutoff, read_qual_cutoff, max_prop_N):
    """
    Parse SAM file contents and return sequence. For matched read pairs,
    merge the reads together into a single sequence
    mapping_cutoff int mapping quality cutoff.  Ignore alignments with mapping quality lower than the cutoff.
    RETURNS:  list of tuples [query name, aligned query seq]
    """
    fasta = []
    with open(sam_filename, 'r') as sam_fh:
        lines = sam_fh.readlines()

        # If this is a completely empty file, return
        if len(lines) == 0:
            return None

        # Skip top SAM header lines
        for start, line in enumerate(lines):
            if not line.startswith('@'):
                break

        # If this is an empty SAM, return
        if start == len(lines) - 1:
            return None

        i = start
        while i < len(lines):
            qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[i].rstrip().split('\t')[:11]

            # If read failed to map or has poor mapping quality, skip it
            if refname == '*' or cigar == '*' or int(mapq) < mapping_cutoff:
                i += 1
                continue

            pos1 = int(pos)
            shift, seq1, qual1 = apply_cigar(cigar, seq, qual)

            if not seq1:
                i += 1
                continue

            seq1 = '-' * pos1 + seq1  # FIXME: We no longer censor bases up front
            qual1 = '-' * pos1 + qual1  # FIXME: Give quality string the same offset


            # No more lines
            if (i + 1) == len(lines):
                break

            # Look ahead in the SAM for matching read
            qname2, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[i + 1].rstrip().split('\t')[:11]

            if qname2 == qname:

                # Second read failed to map
                if refname == '*' or cigar == '*':
                    fasta.append([qname, seq1])
                    i += 2
                    continue

                pos2 = int(pos)
                shift, seq2, qual2 = apply_cigar(cigar, seq, qual)

                # Failed to parse CIGAR
                if not seq2:
                    fasta.append([qname, seq1])
                    i += 2
                    continue

                #seq2 = '-'*pos2 + censor_bases(seq2, qual2, cutoff)
                seq2 = '-' * pos2 + seq2  # FIXME: We no longer censor bases up front
                qual2 = '-' * pos2 + qual2  # FIXME: Give quality string the same offset

                mseq = merge_pairs(seq1, seq2, qual1, qual2, read_qual_cutoff)  # FIXME: We now feed these quality data into merge_pairs

                # Sequence must not have too many censored bases
                if mseq.count('N') / float(len(mseq)) < max_prop_N:
                    fasta.append([qname, mseq])

                i += 2
                continue

            # ELSE no matched pair
            fasta.append([qname, seq1])
            i += 1

    return fasta


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

    return (res)


def create_depth_file_from_bam (bam_filename):
    """
    Gets the coverage from samtools depth.
    Creates a sorted bam file with the same name as bam_filename but appended with ".sort".
    Creates a samtools per-base depth file with the same name as bam_filename but appended with ".depth".
    str bam_filename:  full path to bam file
    Raises subprocess.CalledProcessError
    """

    import subprocess

    # Get per-base depth
    depth_filename = bam_filename + ".depth"
    with open(depth_filename, 'w') as depth_fh:
        # Sort the bam file by leftmost position on the reference assembly.  Required for samtools depth.
        sorted_bam_filename = bam_filename + ".sort"
        subprocess.check_call('samtools', 'sort', '-f', bam_filename, sorted_bam_filename, stderr=subprocess.STDOUT, shell=False)
        depth_output = subprocess.check_output(['samtools', 'depth', bam_filename, '>', ], stderr=subprocess.STDOUT, shell=False)


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
