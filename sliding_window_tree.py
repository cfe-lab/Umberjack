import slice_miseq
import sam_handler
import os
import subprocess
import Utility

# TODO:  do not hardcode these
SELECTION_BF = "../hyphy/getDnDs.bf"

HYPHY_EXE = "/home/thuy/programs/hyphy/hyphy-2.2/HYPHYMP"
HYPHY_LIBDIR = "/home/thuy/programs/hyphy/hyphy-2.2/res/TemplateBatchFiles"
HYPHY_BASEDIR = "/home/thuy/programs/hyphy/hyphy-2.2"

FASTTREE_EXE = "/home/thuy/programs/fasttree/FastTree"


def process_windows(ref_fasta_filename, sam_filename, mapping_cutoff, read_qual_cutoff, max_prop_N,
                    window_depth_thresh, window_breadth_thresh):
    """
    Slides window along genome.
    Creates the multiple sequence aligned fasta files for the window.
    Feeds the multiple-sequence aligned fasta files to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.
    """
    # TODO:  how to get ORFs?
    # TODO:  do not use kalign2 from PATH
    # TODO:  do not hyphy from PATH

    # keep track of reference contig/chromosome and its length
    ref_to_len = Utility.get_seq2len(fasta_filename=ref_fasta_filename)

    # For every contig/chromosome in the reference fasta, create a pseudo multiple-sequence aligned fasta file
    # TODO:  handle indels in multiple sequence align file
    # TODO:  allow option to specify msa file location
    sam_filename_prefix, fileExtension = os.path.splitext(sam_filename)
    msa_fasta_filename = sam_filename_prefix + ".msa.fasta"

    # TODO:  need to create separate msa fasta file for each reference contig/chromosome
    sam_handler.get_msa_fasta_from_sam(sam_filename=sam_filename, ref_fasta_filename=ref_fasta_filename,
                                       mapping_cutoff=mapping_cutoff, read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N,
                                       out_fasta_filename=msa_fasta_filename)

    #  Get the best window size
    # TODO:  handle multiple contig/chromosome in reference fasta.  Right now only handles one contig/chromosome.
    # windowsize = slice_miseq.get_best_window_size_from_sam(sam_filename=sam_filename, ref_filename=ref_fasta_filename,
    #                                             depth_thresh=window_depth_thresh, breadth_thresh=window_breadth_thresh)

    windowsize = 300    # TODO:  do not hardcode this
    for ref in ref_to_len:
        last_window_start_nucpos = ref_to_len[ref] - windowsize
        for start_nucpos in range(1, last_window_start_nucpos):  # start_nucpos is 1-based nucleotide position
            if ((start_nucpos - 1) % 3) != 0:
                continue    # Windows start on the beginning of the codon

            end_nucpos = start_nucpos + windowsize                # end_nucpos is 1-based nucleotide position

            # Slice the multiple sequence aligned fasta file into a window fasta
            msa_window_fasta_filename = slice_miseq.create_slice_msa_fasta(fasta_filename=msa_fasta_filename, start_pos=start_nucpos, end_pos=end_nucpos)

            # Feed window fasta into fasttree to make a tree
            window_prefix = os.path.splitext(msa_window_fasta_filename)[0]
            fastree_logfilename = window_prefix + ".fasttree.log"
            fastree_treefilename = window_prefix + ".tree"
            fasttree_stdouterr_filename = window_prefix + ".fasttree.stdouterr.txt"
            with open(fasttree_stdouterr_filename, 'w') as fasttree_stdouterr_fh:
                subprocess.check_call([FASTTREE_EXE, '-gtr', '-nt', '-nosupport',
                                       '-log', fastree_logfilename, '-out', fastree_treefilename,
                                       msa_window_fasta_filename],
                                      stdout=fasttree_stdouterr_fh, stderr=fasttree_stdouterr_fh, shell=False)


            hyphy_modelfit_filename = window_prefix + ".nucmodelfit"
            hyphy_dnds_tsv_filename = window_prefix + ".dnds.tsv"
            pvalue = 0.95   # TODO:  do not hardcode
            hyphy_input_str = "\n".join(["Universal",
                                        "New Analysis",
                                        os.path.abspath(msa_window_fasta_filename),
                                        "Default",
                                        os.path.abspath(fastree_treefilename),
                                        os.path.abspath(hyphy_modelfit_filename),
                                        "Neutral",
                                        "Single Ancestor Counting",
                                        "-1",
                                        "Full tree",
                                        "Averaged",
                                        "Approximate",
                                        "0.95",
                                        "Export to File",
                                        os.path.abspath(hyphy_dnds_tsv_filename),
                                        "Count\n"])

            # Feed window tree into hyphy to find dnds for the window
            hyphy_log = window_prefix + ".hyphy.log"
            with open(hyphy_log, 'w') as hyphy_out_fh:
                hyphy_proc = subprocess.Popen([HYPHY_EXE, "LIBPATH", HYPHY_LIBDIR, "BASEPATH", HYPHY_BASEDIR, SELECTION_BF],
                                              stdin=subprocess.PIPE, stdout=hyphy_out_fh, stderr=hyphy_out_fh, shell=False)
                hyphy_proc.communicate(hyphy_input_str)

    dnds_tsv_dir = os.path.dirname(hyphy_dnds_tsv_filename)
    ref2SeqDnDsInfo = slice_miseq.get_seq_dnds(dnds_tsv_dir=dnds_tsv_dir, pvalue_thresh=pvalue, ref_fasta_filename=ref_fasta_filename)

    return ref2SeqDnDsInfo



