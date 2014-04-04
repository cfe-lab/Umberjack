import sam_handler
import os
import subprocess
import StringIO


# TODO:  split MSA fasta files by reference contig/chromosome
def slice_msa_fasta(fasta_filename, start_pos, end_pos):
    """
    From a fasta file of multiple sequence alignments, extract the sequence sequence from desired region.
    Writes the extracted sequences into a new fasta file with ".<start_pos>_<end_pos>.fasta" suffix.
    If the sequence is shorter than <end_pos>, then fills in the gaps with '-' characters so that it ends at <end_pos>
    :rtype : object list of tuples [[header1, seq1], [header2, seq2]...]
    :param fasta_filename: full file path to fasta of multiple sequence alignments
    :param start_pos: 1-based start position of region to extract
    :param end_pos: 1-based end position of region to extract
    """
    fasta_filename_prefix, fileExtension = os.path.splitext(fasta_filename)
    slice_fasta_filename = fasta_filename_prefix + "." + str(start_pos) + "_" + str(end_pos) + ".fasta"
    with open(fasta_filename, 'r') as fasta_fh:
        with open(slice_fasta_filename, 'w') as slice_fasta_fh:
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


def get_fasta_headers(fasta_filename):
    """
    Gets a list of headers from the fasta.  Does not include the ">" header prefix.
    Does not include anything after the first whitespace in the header.
    :rtype list[str] : list of headers
    :param fasta_filename : full file path to the fasta file
    """
    headers = []
    with open(fasta_filename, 'r') as fasta_fh:
        for line in fasta_fh:
            if line[0] == '>':
                header = line[1:].rstrip().split()
                headers.append(header)
    return headers




def get_best_window_size(sam_filename, ref_filename, depth_thresh, breadth_thresh):
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
    """

    # convert sam to bam
    bam_filename = sam_handler.sam_to_sort_bam(sam_filename=sam_filename, ref_filename=ref_filename)

    # create depth file for all references
    depth_filename = sam_handler.create_depth_file_from_bam(bam_filename=bam_filename)

    best_size = get_best_window_size_from_depth(depth_filename, depth_thresh, breadth_thresh)
    return best_size


def get_longest_seq_size_from_fasta(fasta_filename):
    """
    Gets the size of the longest sequence in the fasta.
    :rtype int :  size in bp of the longest sequence in the fasta.  Or -1 if error.
    :param fasta_filename: full filepath to the fasta.
    """
    longest_seq_len = -1
    with open(fasta_filename, 'r') as fasta_fh:
        seq_len = 0
        for line in fasta_fh:
            line = line.rstrip()
            if line[0] == '>':
                longest_seq_len = max(seq_len, longest_seq_len)
                seq_len = 0
            else:
                seq_len += len(line)
        longest_seq_len = max(seq_len, longest_seq_len)

    return longest_seq_len


def is_valid_window_exist_from_depth(depth_filename, ref_fasta_filename, depth_thresh, breadth_thresh):
    """
    In order for a window to be valid, it must be have at least depth_thresh reads that cover
    <breadth_cov_thresh> fraction of the window.

    The largest window size possible is the size of the longest reference contig.
    If the longest stretch of consecutive bases with unacceptable depth coverage
    is larger than the breadth_thres * largest possible window size,
    then there are no valid windows.

    :rtype boolean: whether a valid window exists
    :param depth_filename:  full path to samtools depth file
    :param ref_fasta_filename: full path to reference fasta file
    :param depth_thresh: minimum reads that must exist in the window that cover breadth_thresh fraction of the window.
    :param breadth_thresh:  minimum fraction of window that must be covered by at least depth_thresh reads.
    """
    is_valid = False
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

        max_poss_windowsize = get_longest_seq_size_from_fasta(fasta_filename=ref_fasta_filename)

        if longest_consec_base_below_depth <= (breadth_thresh * max_poss_windowsize):
            is_valid = True

    return is_valid


def get_best_window_size_from_depth(sorted_sam_filename, depth_filename, ref_fasta_filename, depth_thresh, breadth_thresh):
    """
    Slide through the alignments for every sequence in the reference file.
    Find the window size such that the smallest window size is selected such that
    for every window of that size along the genome,
    the window is covered by at least <min_reads> that cover <breadth_cov_thresh> fraction of the window.

    :rtype int :  best window size in bases or -1 if there is no window size that meets the constraints
    :param str sorted_sam_filename: full file path to sam file sorted by left coordinates
    :param str depth_filename: full file path to the samtools depth file
    :param str ref_fasta_filename: full file path to the reference fasta file
    :param float breadth_thresh : fraction of window that reads must cover
    :param int depth_thresh: minimum number of reads in window that must cover <breadth_thresh> fraction of the window
    """

    best_size = 1

    # Do a pass to check that there is a window size that meets the constraints
    is_valid = is_valid_window_exist_from_depth(depth_filename=depth_filename, ref_fasta_filename=ref_fasta_filename,
                                     depth_thresh=depth_thresh, breadth_thresh==breadth_thresh)

    if not is_valid:
        return -1

    # get end position of alignment based on cigar
    sam_output = StringIO.StringIO()
    # TODO:  how to pipe line by line to stringio from samtools?????
    subprocess.check_call(['samtools', 'view', sorted_sam_filename], stdout=sam_output, shell=False)




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
