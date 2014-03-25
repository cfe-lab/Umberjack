import sam_handler
import os

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




def get_best_window_size(sam_filename, ref_filename, cov_thresh):
    """
    Create index reference file.
    Create sorted, indexed bam file from sam file.
    Create depth file from bam file.
    Slide through the alignments for every sequence in the reference file.  Find the window size such that the average coverage is
    over the given coverage threshold for all sequences in the reference file.

    :rtype int :  best window size
    :param sam_filename: full file path to the sam alignment file
    :param ref_filename: full file path to the reference file
    :param cov_thresh: the minimum average read coverage in window
    """

    # convert sam to bam
    bam_filename = sam_handler.sam_to_sort_bam(sam_filename=sam_filename, ref_filename=ref_filename)

    # create depth file for all references
    depth_filename = sam_handler.create_depth_file_from_bam(bam_filename=bam_filename)

    best_size =  get_best_window_size_from_depth(depth_filename, cov_thresh)
    return best_size


def get_best_window_size_from_depth(depth_filename, cov_thresh):
    """
    Slide through the alignments for every sequence in the reference file.  Find the window size such that the average coverage is
    over the given coverage threshold for all sequences in the reference file.

    :rtype int :  best window size
    :param depth_filename: full file path to the samtools depth file
    :param cov_thresh: the minimum average read coverage in window
    """

    best_size = 1
    # iterate through each reference sequence and find the best window size across all ref seq
    with open(depth_filename, 'r') as depth_fh:
        # assume that depth file is sorted by reference, then position
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
                if ave_window_cov < cov_thresh:
                    raise Exception("Reference " + lastref + " average coverage of " + ave_window_cov +
                                    " does not meet the minimum coverage threshold of " + cov_thresh)
                windowsize = 0
                total_window_cov = 0
                last_window_end = 0
                lastref = ""

            # continue sliding along current ref
            windowsize = int(pos) - last_window_end
            total_window_cov += int(depth)
            ave_window_cov = float(total_window_cov) / windowsize

            if ave_window_cov < cov_thresh:  # extend window
                best_size = max(best_size, windowsize)
            else:  # window has sufficient coverage.  Start another window.
                windowsize = 0
                total_window_cov = 0
                last_window_end = int(pos)

            lastref = ref

        return best_size  # NOPE!  check out http://seqanswers.com/forums/showthread.php?t=25587 ?