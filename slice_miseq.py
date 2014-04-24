import sam_handler
import os
import subprocess
import StringIO
import Utility
import csv
import glob


HYPHY_TSV_SCALED_DNDS_COL = 'Scaled dN-dS'
HYPHY_TSV_PROB_FULLSEQ_NS = 'P{S leq. observed}'

NUC_PER_CODON = 3


class SiteDnDsInfo:
    def __init__(self):
        self.accum_win_dnds = 0
        self.total_win_cover_site = 0

    def add_dnds(self, dnds):
        self.total_win_cover_site += 1
        self.accum_win_dnds += dnds

    def get_ave_dnds(self):
        if self.accum_win_dnds == 0 and self.total_win_cover_site == 0:
            return None
        else:
            return self.accum_win_dnds / self.total_win_cover_site


class SeqDnDsInfo:
    """
    A sequence of Dn/Ds by site
    """
    def __init__(self, seq_len):
        self.dnds_seq = [SiteDnDsInfo()] * seq_len

    def add_site_dnds(self, site, dnds):
        self.dnds_seq[site].add_dnds(dnds)

    def get_site_ave_dnds(self, site):
        return self.dnds_seq[site].get_ave_dnds()

    def get_seq_len(self):
        return len(self.dnds_seq)







# TODO:  finish me!!!
def create_slice_msa_fasta_from_bam(bam_filename, start_pos, end_pos):
    """
    From a fasta file of multiple sequence alignments, extract the sequence sequence from desired region.
    Writes the extracted sequences into a new fasta file with ".<start_pos>_<end_pos>.fasta" suffix.
    If the sequence is shorter than <end_pos>, then fills in the gaps with '-' characters so that it ends at <end_pos>

    ASSUMES:  that the start and end position of the clipped, unpadded sequence is in the header with this format
    >sequence_name start_pos_1based end_pos_1based
    :rtype str: full filepath to sliced multiple sequence alignment fasta
    :param bam_filename: full file path to bam file sorted by left coordinate
    :param start_pos: 1-based start position of region to extract
    :param end_pos: 1-based end position of region to extract
    """
    fasta_filename_prefix, fileExtension = os.path.splitext(fasta_filename)
    slice_fasta_filename = fasta_filename_prefix + "." + str(start_pos) + "_" + str(end_pos) + ".fasta"
    with open(fasta_filename, 'r') as fasta_fh, open(slice_fasta_filename, 'w') as slice_fasta_fh:
        header = ""
        seq = ""
        for line in fasta_fh:
            line = line.rstrip()  # remove trailing whitespace

            if line[0] == '>':  # last sequence is finished.  Write out last sequence
                if seq:
                    pad = ''
                    if end_pos > len(seq):
                        pad = '-' * (end_pos - len(seq))
                    slice_fasta_fh.write(seq[start_pos-1:end_pos] + pad + "\n")


                seq = ""

                header = line   # Write out current header
                slice_fasta_fh.write(header + "\n")
            else:   # cache current sequence so that entire sequence is on one line
                seq += line

        if seq:   # end of fasta file, write out the last sequence still in cache
            pad = ''
            if end_pos > len(seq):
                pad = '-' * (end_pos - len(seq))
            slice_fasta_fh.write(seq[start_pos-1:end_pos] + pad + "\n")
    return slice_fasta_filename


def __get_window_seq(seq, start_pos, end_pos, breadth_thresh=0.0):
    """
    If the sequence fits the window constraints, then returns the sliced sequence.
    Otherwise returns None.

    :rtype str :  Sequence that fits in the window, right padded if necessary.  None if sequence doesn't fit window constraints.
    :param seq str: sequence with no newlines
    :param start_pos int: 1-based
    :param end_pos int : 1-based
    :param breadth_thresh float: fraction of sequence that be A,C,T,or G within start_pos and end_pos inclusive.
    """
    if seq and len(seq) >= start_pos:
        slice_seq = seq[start_pos-1:end_pos]
        slice_seq_uc = slice_seq.upper()
        a_count = slice_seq_uc.count('A')
        t_count = slice_seq_uc.count('T')
        g_count = slice_seq_uc.count('G')
        c_count = slice_seq_uc.count('C')
        if float(a_count+g_count+c_count+t_count)/(end_pos - start_pos + 1) >= breadth_thresh:
            pad = '-' * max(0, end_pos - len(seq))
            return slice_seq

    return None


def create_slice_msa_fasta(fasta_filename, out_fasta_filename, start_pos, end_pos, breadth_thresh=0.0):
    """
    From a fasta file of multiple sequence alignments, extract the sequence sequence from desired region.
    Writes the extracted sequences into a new fasta file with ".<start_pos>_<end_pos>.fasta" suffix.
    If the sequence is shorter than <end_pos>, then fills in the gaps with '-' characters so that it ends at <end_pos>

    Only puts in the read into the sliced msa fasta if it obeys the window constraints.

    :rtype int: total sequences written
    :param str fasta_filename: full file path to fasta of multiple sequence alignments
    :param str out_fasta_filename:  full file path to output fasta of sliced multiple sequence alignments
    :param int start_pos : 1-based start position of region to extract
    :param int end_pos: 1-based end position of region to extract
    :param float breadth_thresh: fraction of sequence that be A,C,T,or G within start_pos and end_pos inclusive.
    """

    total_seq = 0
    with open(fasta_filename, 'r') as fasta_fh, open(out_fasta_filename, 'w') as slice_fasta_fh:
        header = ""
        seq = ""
        for line in fasta_fh:
            line = line.rstrip().split()[0]  # remove trailing whitespace and any test after the first whitespace

            if line[0] == '>':  # previous sequence is finished.  Write out previous sequence
                window_seq = __get_window_seq(seq=seq, start_pos=start_pos, end_pos=end_pos, breadth_thresh=breadth_thresh)
                if window_seq:
                    slice_fasta_fh.write(header + "\n")
                    slice_fasta_fh.write(window_seq + "\n")
                    total_seq += 1

                seq = ""
                header = line

            else:   # cache current sequence so that entire sequence is on one line
                seq += line

        window_seq = __get_window_seq(seq=seq, start_pos=start_pos, end_pos=end_pos, breadth_thresh=breadth_thresh)
        if window_seq:
            slice_fasta_fh.write(header + "\n")
            slice_fasta_fh.write(window_seq + "\n")
            total_seq += 1

    return total_seq


def get_best_window_size_from_sam(sam_filename, ref_filename, depth_thresh, breadth_thresh, min_size):
    """
    Create index reference file.
    Create sorted, indexed bam file from sam file.
    Create depth file from bam file.
    Slide through the alignments for every sequence in the reference file.
    The smallest window size is selected such that for every window of that size along the genome,
    the window is covered by at least <min_reads> that cover <breadth_cov_thresh> fraction of the window.

    :rtype int :  best window size in bases or -1 if there are none that meet the constraints
    :param sam_filename: full file path to the sam alignment file
    :param ref_filename: full file path to the reference file
    :param float breadth_thresh : fraction of window that reads must cover
    :param int depth_thresh: minimum number of reads in window that must cover <breadth_cov_thresh> fraction of the window
    :param int min_size: minimum number of bases allowed in the window.
    """

    # convert sam to bam
    bam_filename = sam_handler.sam_to_sort_bam(sam_filename=sam_filename, ref_filename=ref_filename)

    # create depth file for all references
    depth_filename = sam_handler.create_depth_file_from_bam(bam_filename=bam_filename)

    best_size = get_best_window_size_from_depthfile(sorted_bam_filename=bam_filename,
                                                    depth_filename=depth_filename,
                                                    ref_fasta_filename=ref_filename,
                                                    depth_thresh=depth_thresh, breadth_thresh=breadth_thresh, min_size=min_size)

    return best_size


def get_longest_consec_base_below_depth(depth_filename, ref_fasta_filename, depth_thresh):
    """
    In order for a window to be valid, it must be have at least depth_thresh reads that cover
    <breadth_cov_thresh> fraction of the window.

    The largest window size possible is the size of the longest reference contig.
    If the longest stretch of consecutive bases with unacceptable depth coverage
    is larger than the breadth_thres * largest possible window size,
    then there are no valid windows.

    :rtype int: total consecutive bases below the specified depth
    :param depth_filename:  full path to samtools depth file
    :param ref_fasta_filename: full path to reference fasta file
    :param depth_thresh: minimum reads that must exist in the window that cover breadth_thresh fraction of the window.
    :param breadth_thresh:  minimum fraction of window that must be covered by at least depth_thresh reads.
    """
    # is_valid = False
    with open(depth_filename, 'r') as depth_fh:

        # Assume that depth file is sorted by reference, then position
        # Assume position is 1-based.
        # Assume that the depth file skips bases with zero coverage.
        last_ref = ''
        last_pos = 0
        longest_consec_base_below_depth = 0
        consec_base_below_depth = 0
        for line in depth_fh:
            ref, pos, depth = line.rstrip().split()  # pos is 1-based
            if depth < depth_thresh:
                if last_ref != ref:
                    longest_consec_base_below_depth = max(longest_consec_base_below_depth, consec_base_below_depth)
                    consec_base_below_depth = 0

                consec_base_below_depth += (pos - last_pos)

            else:
                longest_consec_base_below_depth = max(longest_consec_base_below_depth, consec_base_below_depth)
                consec_base_below_depth = 0

            last_ref = ref
            last_pos = pos

        max_poss_windowsize = Utility.get_longest_seq_size_from_fasta(fasta_filename=ref_fasta_filename)

        # if longest_consec_base_below_depth <= (breadth_thresh * max_poss_windowsize):
        #     is_valid = True

    return longest_consec_base_below_depth


def is_window_size_valid(depth_filename, ref_fasta_filename, depth_thresh, breadth_thresh, window_size):

    with open(depth_filename, 'r') as depth_fh:

        # Assume that depth file is sorted by reference, then position
        # Assume position is 1-based.
        # Assume that the depth file skips bases with zero coverage.
        last_ref = ''
        last_pos = 0
        longest_consec_base_below_depth = 0
        consec_base_below_depth = 0
        for line in depth_fh:
            ref, pos, depth = line.rstrip().split()  # pos is 1-based
            if depth < depth_thresh:
                if last_ref != ref:
                    longest_consec_base_below_depth = max(longest_consec_base_below_depth, consec_base_below_depth)
                    consec_base_below_depth = 0

                consec_base_below_depth += (pos - last_pos)

            else:
                longest_consec_base_below_depth = max(longest_consec_base_below_depth, consec_base_below_depth)
                consec_base_below_depth = 0

            last_ref = ref
            last_pos = pos

def get_best_window_size_from_depthfile(sorted_bam_filename, depth_filename, ref_fasta_filename, depth_thresh, breadth_thresh, min_size):
    """
    Slide through the alignments for every sequence in the reference file.
    Find the window size such that the smallest window size is selected such that
    for every window of that size along the genome,
    the window is covered by at least <min_reads> that cover <breadth_cov_thresh> fraction of the window.

    :rtype int :  best window size in bases or -1 if there is no window size that meets the constraints
    :param str sorted_bam_filename: full file path to bam file sorted by left coordinates
    :param str depth_filename: full file path to the samtools depth file
    :param str ref_fasta_filename: full file path to the reference fasta file
    :param float breadth_thresh : fraction of window that reads must cover
    :param int depth_thresh: minimum number of reads in window that cover at least <breadth_thresh> fraction of the window
    :param int min_size: minimum number of bases allowed in the window.
    """

    best_size = 1

    # Find the smallest possible window size
    size_bad_cov = get_longest_consec_base_below_depth(depth_filename=depth_filename, ref_fasta_filename=ref_fasta_filename,
                                                   depth_thresh=depth_thresh)
    windowsize = size_bad_cov/(1-breadth_thresh)  # our initial window size.

    # get end position of alignment based on cigar
    # sam_output = StringIO.StringIO()
    # # TODO:  how to pipe line by line to stringio from samtools?????
    # subprocess.check_call(['samtools', 'view', sorted_bam_filename], stdout=sam_output, shell=False)



    # keep track of positions with indels
    with open(depth_filename, 'r') as depth_fh:
        lastref = ""
        last_window_end = 0
        windowsize = 0
        total_window_cov = 0
        for line in depth_fh:
            ref, pos, depth = line.rstrip().split()  # pos is 1-based

            # we slid past the last ref
            if not ref == lastref and lastref:
                ave_window_cov = float(total_window_cov) / (windowsize + 1)

                # insufficient coverage over the entire last ref
                if ave_window_cov < depth_thresh:
                    raise Exception("Reference " + lastref + " average coverage of " + ave_window_cov +
                                    " does not meet the minimum coverage threshold of " + depth_thresh)
                windowsize = 0
                total_window_cov = 0
                last_window_end = 0
                lastref = ""

            # continue sliding along current ref
            windowsize = int(pos) - last_window_end
            total_window_cov += int(depth)
            ave_window_cov = float(total_window_cov) / windowsize

            if ave_window_cov < depth_thresh:  # extend window
                best_size = max(best_size, windowsize)
            else:  # window has sufficient coverage.  Start another window.
                windowsize = 0
                total_window_cov = 0
                last_window_end = int(pos)

            lastref = ref

        return best_size  # NOPE!  check out http://seqanswers.com/forums/showthread.php?t=25587 ?


def get_seq_dnds(dnds_tsv_dir, ref_fasta_filename, pvalue_thresh):
    """
    For an entire window, get the average dn/ds for every base position.
    We expect that all HyPhy has written out dn/ds tsv files for every window.

    :rtype dict of SeqDnDsInfo
    :param dnds_tsv_dir: 
    :param ref_fasta_filename: 
    :param pvalue_thresh:

    Assumes these are the columns in the HyPhy DN/DS TSV output:
    Observed S Changes
    Observed NS Changes
    E[S Sites]: proportion of random one-nucleotide substitutions that are expected to be synonymous
    E[NS Sites]: roportion of random one-nucleotide substitutions that are expected to be non-synonymous
    Observed S. Prop.
    P{S}:  proportion of substitutions expected to be synonymous under neutral evolution
    dS
    dN
    dN-dS
    P{S leq. observed}:  binomial distro pvalue.  Probability of getting less than the observed synynomous substitutions
        under the binomial distribution where probability of 1 synonymous codon = P{S}
    P{S geq. observed}: 1-pvalue.  Probability of getting more than the observed synynomous substitutions
        under the binomial distribution where probability of 1 synonymous codon = P{S}
    Scaled dN-dS:  dN-dS normalized by the total length of the tree.

    We only want significant Dn/Ds.  Use the normalized value so that dn/ds can be compared from every window.
    """

    refs = Utility.get_fasta_headers(fasta_filename=ref_fasta_filename)
    ref2nuclen = Utility.get_seq2len(ref_fasta_filename)
    ref2SeqDnDsInfo = {}
    for dnds_tsv_filename in glob.glob(dnds_tsv_dir + os.sep + "*.dnds.tsv"):
        with open(dnds_tsv_filename, 'r') as dnds_fh:

            for ref in refs:
                if not ref in ref2SeqDnDsInfo:
                    # TODO:  check what happens if window length is not divisible by 3
                    ref_codon_len = ref2nuclen[ref]/NUC_PER_CODON
                    ref2SeqDnDsInfo[ref] = SeqDnDsInfo(ref_codon_len)

                # TestSample-RT_S17.HIV1B-vif.remap.msa.14_114.dnds.tsv
                # NB:  the *dnds.tsv file names use 1-based nucleotide position numbering
                # Window starts at this 1-based nucleotide position with respect to the reference
                win_start_nuc_pos_1based_wrt_ref = int(dnds_tsv_filename.split('.dnds.tsv')[0].split('.')[-1].split('_')[0])
                # Window starts at this 0-based codon position with respect to the reference
                win_start_codon_idx_0based_wrt_ref = (win_start_nuc_pos_1based_wrt_ref - 1)/NUC_PER_CODON

                reader = csv.DictReader(dnds_fh, delimiter='\t',)
                for offset_0based, codon_row in enumerate(reader):    # Every codon site is a row in the *.dnds.tsv file
                    pval = float(codon_row[HYPHY_TSV_PROB_FULLSEQ_NS])
                    if pval <= pvalue_thresh:
                        ref_codon_idx_0based = win_start_codon_idx_0based_wrt_ref + offset_0based

                        # Dict  {ref1:SeqDnDsInfo, ref2: SeqDnDsInfo}
                        dnds = float(codon_row[HYPHY_TSV_SCALED_DNDS_COL])
                        ref2SeqDnDsInfo[ref].add_site_dnds(ref_codon_idx_0based, dnds)

        return ref2SeqDnDsInfo



