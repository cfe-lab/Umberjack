import slice_miseq
import sam_handler
import os, sys
import subprocess
import Utility
import logging
import errno

# TODO:  do not hardcode these
# SELECTION_BF = "../hyphy/getDnDs.bf"
SELECTION_BF = "res/TemplateBatchFiles/QuickSelectionDetection.bf"
HYPHY_EXE = "/home/thuy/programs/hyphy/hyphy-2.2/HYPHYMP"
HYPHY_LIBDIR = "/home/thuy/programs/hyphy/hyphy-2.2/res/TemplateBatchFiles/"
HYPHY_BASEDIR = "/home/thuy/programs/hyphy/hyphy-2.2/"

FASTTREE_EXE = "/home/thuy/programs/fasttree/FastTree"

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)

NUC_PER_CODON = 3
        

def process_windows(ref_fasta_filename, sam_filename, out_dir, mapping_cutoff, read_qual_cutoff, max_prop_N,
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
    ref_to_nuclen = Utility.get_fasta_headers(fasta_filename=ref_fasta_filename)

    # For every reference contig/chromosome, find the ave dnds by codon site
    for ref in ref_to_nuclen:
        ref_out_dir = out_dir + os.sep + ref
        Utility.create_dir_check(ref_out_dir)

        # Create a pseudo multiple-sequence aligned fasta file
        # TODO:  handle indels in multiple sequence align file
        sam_filename_nopath = os.path.split(sam_filename)[1]
        sam_filename_prefix = os.path.splitext(sam_filename_nopath)[0]
        msa_filename_prefix = ref_out_dir + os.sep + sam_filename_prefix + "." + ref + ".msa"
        msa_fasta_filename = msa_filename_prefix + ".fasta"
        LOGGER.debug("Start Full MSA-Fasta from SAM for ref " + ref)
        sam_handler.create_msa_fasta_from_sam(sam_filename=sam_filename, ref=ref, ref_len=ref_to_nuclen[ref],
                                              out_fasta_filename=msa_fasta_filename, mapping_cutoff=mapping_cutoff,
                                              read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N)
        LOGGER.debug("Done Full MSA-Fasta from SAM for ref " + ref)


        #  Get the best window size
        # TODO:  handle multiple contig/chromosome in reference fasta.  Right now only handles one contig/chromosome.
        # windowsize = slice_miseq.get_best_window_size_from_sam(sam_filename=sam_filename, ref_filename=ref_fasta_filename,
        #                                             depth_thresh=window_depth_thresh, breadth_thresh=window_breadth_thresh)

        windowsize = 300    # TODO:  do not hardcode this

        last_window_start_nucpos = ref_to_nuclen[ref] - windowsize
        for start_nucpos in range(1, last_window_start_nucpos):  # start_nucpos is 1-based nucleotide position
            if ((start_nucpos - 1) % 3) != 0:
                continue    # Windows start on the beginning of the codon

            end_nucpos = start_nucpos + windowsize                # end_nucpos is 1-based nucleotide position

            # Slice the multiple sequence aligned fasta file into a window fasta
            LOGGER.debug("Start Create Sliced MSA-Fasta " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref)
            msa_window_filename_prefix = msa_filename_prefix + "." + str(start_nucpos) + "_" + str(end_nucpos)
            msa_window_fasta_filename = msa_window_filename_prefix + ".fasta"
            total_slice_seq = slice_miseq.create_slice_msa_fasta(fasta_filename=msa_fasta_filename,
                                                                 out_fasta_filename=msa_window_fasta_filename,
                                                                 start_pos=start_nucpos, end_pos=end_nucpos,
                                                                 breadth_thresh=window_breadth_thresh)
            LOGGER.debug("Done Create Sliced MSA-Fasta " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref +
                         ".  Wrote " + str(total_slice_seq) + " to file")

            if total_slice_seq < window_depth_thresh:
                LOGGER.warn("Window " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref +
                            " does not satisfy window depth constraints")
            else:
                # Feed window fasta into fasttree to make a tree
                LOGGER.debug("Window " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref +
                             " satisfies window depth constraints")
                LOGGER.debug("Start Fasttree for window " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref)
                fastree_logfilename = msa_window_filename_prefix + ".fasttree.log"
                fastree_treefilename = msa_window_filename_prefix + ".tree"
                fasttree_stdouterr_filename = msa_window_filename_prefix + ".fasttree.stdouterr.txt"
                with open(fasttree_stdouterr_filename, 'w') as fasttree_stdouterr_fh:
                    subprocess.check_call([FASTTREE_EXE, '-gtr', '-nt', '-nosupport',
                                           '-log', fastree_logfilename, '-out', fastree_treefilename,
                                           msa_window_fasta_filename],
                                          stdout=fasttree_stdouterr_fh, stderr=fasttree_stdouterr_fh, shell=False)
                LOGGER.debug("Done Fasttree for window " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref)

                LOGGER.debug("Start HyPhy for window " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref)
                hyphy_modelfit_filename = msa_window_filename_prefix + ".nucmodelfit"
                hyphy_dnds_tsv_filename = msa_window_filename_prefix + ".dnds.tsv"
                pvalue = 0.95   # TODO:  do not hardcode

                hyphy_input_str = "\n".join(["1",   # Universal
                                            "1",    # New analysis
                                            os.path.abspath(msa_window_fasta_filename), # codon fasta
                                            "1",    # Substitution Model - Use HK85 and MG94xHKY85
                                            os.path.abspath(fastree_treefilename),      # tree file
                                            os.path.abspath(hyphy_modelfit_filename),   # model fit output file
                                            "1",    # Neutral dN/dS = 1
                                            "1",    # Single Acnestor Counting
                                            "-1",   # Approximate
                                            "1",    # Full tree
                                            "1",    # Averaged
                                            "1",    # Approximate extended binomial distro
                                            "0.95", # pvalue threshold
                                            "2",    # Export to file
                                            os.path.abspath(hyphy_dnds_tsv_filename),   # dN/dS tsv output file
                                            "2\n"]) # Count approximate numbers of dN, dS rate classes supported by data

                # Feed window tree into hyphy to find dnds for the window
                hyphy_log = msa_window_filename_prefix + ".hyphy.log"
                with open(hyphy_log, 'w') as hyphy_log_fh:
                    hyphy_cmd = [HYPHY_EXE, "BASEPATH="+HYPHY_BASEDIR, SELECTION_BF]
                    hyphy_proc = subprocess.Popen(hyphy_cmd,
                                                  stdin=subprocess.PIPE, stdout=hyphy_log_fh, stderr=hyphy_log_fh, shell=False)
                    hyphy_proc.communicate(hyphy_input_str)

                    if hyphy_proc.returncode:
                        raise subprocess.CalledProcessError(cmd=hyphy_cmd, returncode=hyphy_proc.returncode)
                LOGGER.debug("Done HyPhy for window " + str(start_nucpos) + "-" + str(end_nucpos) + " for ref " + ref)

        LOGGER.debug("Start Ave Dn/DS for all windows for ref " + ref)
        dnds_tsv_dir = os.path.dirname(hyphy_dnds_tsv_filename)
        ref2SeqDnDsInfo = slice_miseq.get_seq_dnds_info(dnds_tsv_dir=dnds_tsv_dir, pvalue_thresh=pvalue, ref=ref,
                                                        ref_codon_len=ref_to_nuclen[ref]/NUC_PER_CODON)
        LOGGER.debug("Done Ave Dn/DS for all windows  for ref " + ref)

    return ref2SeqDnDsInfo



