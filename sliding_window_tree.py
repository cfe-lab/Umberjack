import slice_miseq
import sam_handler
import os, sys
import subprocess
import Utility
import logging
import argparse
import pool_traceback
import re
import glob

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)

SELECTION_BF = "QuickSelectionDetection.bf"
NUC_MODEL_CMP_BS = "NucModelCompare.bf"
HYPHY_EXE = "HYPHYMP"
HYPHY_BASEDIR = "/usr/local/lib/hyphy/TemplateBatchFiles/"

FASTTREE_EXE = "FastTreeMP"
ENV_OMP_NUM_THREADS = 'OMP_NUM_THREADS'


# TODO:  handle inserts in multiple sequence align file
# TODO:  handle reads that extend the reference in MSA file
# TODO: this is a performance killing step.  Instead of writing all MSA aligned reads to file, we need to just do it for
#   parts of the genome that are used in the sliding window or hold all of this in memory.
def create_full_msa_fasta(sam_filename, out_dir, ref, ref_len, mapping_cutoff, read_qual_cutoff, max_prop_N):
    """
    Creates a pseudo multiple-sequence aligned fasta file for all reads using pairwise alignment from a SAM file.


    :return: path to multiple sequence aligned fasta file of all reads
    :rtype : str
    :param str sam_filename: filepath to sam file of read alignments to the reference
    :param str out_dir: output directory
    :param str ref: reference name
    :param int ref_len: length of reference in nucleotides
    :param float mapping_cutoff: mapping quality cutoff
    :param float read_qual_cutoff: read quality cutoff.  Bases below this cutoff are converted to N's.
    :param float max_prop_N:  maximum fraction of bases in a merged mate-pair read that is allowed to be N's.
                                Reads that exceed this threshold are thrown out.
    """

    sam_filename_nopath = os.path.split(sam_filename)[1]
    sam_filename_prefix = os.path.splitext(sam_filename_nopath)[0]
    msa_filename_prefix = out_dir + os.sep + sam_filename_prefix + "." + ref + ".msa"
    msa_fasta_filename = msa_filename_prefix + ".fasta"

    LOGGER.debug("Start Full MSA-Fasta from SAM for ref " + ref)
    if not os.path.exists(msa_fasta_filename) or os.path.getsize(msa_fasta_filename) <= 0:
        sam_handler.create_msa_fasta_from_sam(sam_filename=sam_filename, ref=ref, ref_len=ref_len,
                                              out_fasta_filename=msa_fasta_filename, mapping_cutoff=mapping_cutoff,
                                              read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N)
        LOGGER.debug("Done Full MSA-Fasta from SAM for ref " + ref)
    else:
        LOGGER.warn("Found existing Full MSA-Fasta from SAM for ref " + ref + ".  Not regenerating")



    return msa_fasta_filename


def do_hyphy(hyphy_exe, hyphy_basedir, threads, hyphy_filename_prefix, mode, codon_fasta_filename, tree_filename, pvalue):
    if mode == "DNDS":
        hyphy_modelfit_filename = hyphy_filename_prefix + ".nucmodelfit"
        hyphy_dnds_tsv_filename = hyphy_filename_prefix + ".dnds.tsv"
        LOGGER.debug("Start HyPhy for window " + hyphy_dnds_tsv_filename)
        if not os.path.exists(hyphy_dnds_tsv_filename) or os.path.getsize(hyphy_dnds_tsv_filename) <= 0:
            hyphy_input_str = "\n".join(["1",   # Universal
                                        "1",    # New analysis
                                        os.path.abspath(codon_fasta_filename), # codon fasta
                                        "2",    #(2):[Custom] Use any reversible nucleotide model crossed with MG94.
                                        "012345", # GTR
                                        os.path.abspath(tree_filename),      # tree file
                                        os.path.abspath(hyphy_modelfit_filename),   # model fit output file
                                        "3",    #(3):[Estimate] Estimate from data with branch corrections(slower).
                                        "1",    # Single Acnestor Counting
                                        "1",    # Full tree
                                        "1",    # Averaged
                                        "1",    # Approximate extended binomial distro
                                        str(pvalue), # pvalue threshold
                                        "2",    # Export to file
                                        os.path.abspath(hyphy_dnds_tsv_filename),   # dN/dS tsv output file
                                        "2\n"]) # Count approximate numbers of dN, dS rate classes supported by data

            # Feed window tree into hyphy to find dnds for the window
            hyphy_log = hyphy_filename_prefix + ".hyphy.log"
            with open(hyphy_log, 'w') as hyphy_log_fh:
                hyphy_cmd = [hyphy_exe, "BASEPATH="+hyphy_basedir, "CPU=" + str(threads), SELECTION_BF]
                hyphy_proc = subprocess.Popen(hyphy_cmd, stdin=subprocess.PIPE, stdout=hyphy_log_fh, stderr=hyphy_log_fh,
                                              shell=False, env=os.environ)
                hyphy_proc.communicate(hyphy_input_str)

                if hyphy_proc.returncode:
                    raise subprocess.CalledProcessError(cmd=hyphy_cmd, returncode=hyphy_proc.returncode)
            LOGGER.debug("Done HyPhy for window " + hyphy_dnds_tsv_filename)
        else:
            LOGGER.debug("Found existing HyPhy for window " + hyphy_dnds_tsv_filename + ". Not regenerating")
    elif mode == "NUC_SUBST":
        hyphy_modelfit_filename = hyphy_filename_prefix + ".nucmodelfit"
        LOGGER.debug("Start HyPhy for window " + hyphy_modelfit_filename)
        if not os.path.exists(hyphy_modelfit_filename) or os.path.getsize(hyphy_modelfit_filename) <= 0:
            hyphy_input_str = "\n".join([
                                        "2",    # [Global] Model parameters are shared by all branches, branch lengths are estimated independently.
                                        "2",    # [Once] Branch lenghts obtained from the general reversible model are reused for subsequent models.
                                        os.path.abspath(codon_fasta_filename), # codon fasta
                                        os.path.abspath(tree_filename),      # tree file
                                        str(pvalue), # pvalue threshold
                                        "1",    # [No] Do not Save each of the 203 files to a separate file
                                        os.path.abspath(hyphy_modelfit_filename),   # model fit output file
                                        "\n"])

            # Feed window tree into hyphy to find dnds for the window
            hyphy_log = hyphy_filename_prefix + ".hyphy.log"
            with open(hyphy_log, 'w') as hyphy_log_fh:
                hyphy_cmd = [hyphy_exe, "BASEPATH="+hyphy_basedir, "CPU=" + str(threads), NUC_MODEL_CMP_BS]
                hyphy_proc = subprocess.Popen(hyphy_cmd, stdin=subprocess.PIPE, stdout=hyphy_log_fh, stderr=hyphy_log_fh,
                                              shell=False, env=os.environ)
                hyphy_proc.communicate(hyphy_input_str)

                if hyphy_proc.returncode:
                    raise subprocess.CalledProcessError(cmd=hyphy_cmd, returncode=hyphy_proc.returncode)
            LOGGER.debug("Done HyPhy for window " + hyphy_modelfit_filename)
        else:
            LOGGER.debug("Found existing HyPhy for window " + hyphy_modelfit_filename + ". Not regenerating")

        # Hyphy creates a *.nucmodelfit file that contains the best fit model (according to AIC) with this entry.  Parse it.
        #   Model averaged rates relative to AG (REV estimates):
        #       AC =   0.1902	(  0.1781)
        #       AT =   0.2058	(  0.2198)
        #       CG =   0.0573	(  0.0567)
        #       CT =   1.2453	(  1.2953)
        #       GT =   0.4195	(  0.4246)
        with open(hyphy_modelfit_filename, 'rU') as fh_nucmodelfit, open(hyphy_modelfit_filename + ".csv", 'w') as fh_nucmodelcsv:
            is_found_rates = False
            fh_nucmodelcsv.write("StartBase,EndBase,Rate\n")
            for line in fh_nucmodelfit:
                if not is_found_rates and "Model averaged rates relative to AG (REV estimates)" in line:
                    is_found_rates = True
                    continue

                if is_found_rates:
                    line = line.rstrip().lstrip()
                    match = re.findall(r'([A-Z][A-Z])\s*=\s*(\d+\.\d+)\s*\(\s*(\d+\.\d+)\s*\)', line, re.IGNORECASE)
                    if not match:
                        break
                    subst, rate, reverse_rate = match[0] # list of 1 tuple
                    init_base, end_base = list(subst)
                    fh_nucmodelcsv.write(init_base + "," + end_base + "," + rate + "\n")
                    fh_nucmodelcsv.write(end_base + "," + init_base + "," + reverse_rate + "\n")
    else:
        raise ValueError("Invalid mode=" + mode)


# TODO:  do multiple test corrections for pvalues
def eval_window(msa_fasta_filename, window_depth_thresh, window_breadth_thresh, start_nucpos, end_nucpos, pvalue, threads, mode="DNDS",
                hyphy_exe=HYPHY_EXE, hyphy_basedir=HYPHY_BASEDIR, fastree_exe=FASTTREE_EXE):
    """
    Handles the processing for a single window along the genome.
    Creates the multiple sequence aligned fasta file for the window.
    Feeds the window multiple-sequence aligned fasta file to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.

    :param str msa_fasta_filename: full filepath to multiple sequence aligned file for all reads.
    :param int window_depth_thresh:  the minimum number of required reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_thresh: the minimum fraction of a window that merged paired-end read must cover to be included in the window.
    :param int start_nucpos:  1-based start nucleotide position of the window
    :param int end_nucpos:  1-based end nucleotide position of the window
    :param float pvalue:  pvalue threshold for detecting dN/dS (used by HyPhy)
    :param int threads: number of threads allotted to processing this window  (only FastTree and HyPhy will be multithreaded)
    :param str hyphy_exe: full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe: full filepath to FastTreeMP executable
    """

    LOGGER.debug("msa_fasta_filename=" + msa_fasta_filename + "\n" +
                 "window_depth_thresh=" + str(window_depth_thresh) + "\n" +
                 "window_breadth_thresh=" + str(window_breadth_thresh) + "\n" +
                 "start_nucpos=" + str(start_nucpos) + "\n" +
                 "end_nucpos=" + str(end_nucpos) + "\n" +
                 "pvalue=" + str(pvalue) + "\n" +
                 "threads=" + str(threads) + "\n")

    # Slice the multiple sequence aligned fasta file into a window fasta
    msa_fasta_filename_prefix = os.path.splitext(msa_fasta_filename)[0]
    msa_window_filename_prefix = msa_fasta_filename_prefix + "." + str(start_nucpos) + "_" + str(end_nucpos)
    msa_window_fasta_filename = msa_window_filename_prefix + ".fasta"
    total_slice_seq = -1
    LOGGER.debug("Start Create Sliced MSA-Fasta " + msa_window_fasta_filename)
    if not os.path.exists(msa_window_fasta_filename) or os.path.getsize(msa_window_fasta_filename) <= 0:
        total_slice_seq = slice_miseq.create_slice_msa_fasta(fasta_filename=msa_fasta_filename,
                                                             out_fasta_filename=msa_window_fasta_filename,
                                                             start_pos=start_nucpos, end_pos=end_nucpos,
                                                             breadth_thresh=window_breadth_thresh)
        LOGGER.debug("Done Create Sliced MSA-Fasta " + msa_window_fasta_filename +
                     ".  Wrote " + str(total_slice_seq) + " to file")
    else:
        LOGGER.warn("Found existing Sliced MSA-Fasta " + msa_window_fasta_filename + ". Not regenerating.")

    fastree_logfilename = msa_window_filename_prefix + ".fasttree.log"
    fastree_treefilename = msa_window_filename_prefix + ".tree"
    fasttree_stdouterr_filename = msa_window_filename_prefix + ".fasttree.stdouterr.txt"
    LOGGER.debug("Start Fasttree for window " + fastree_treefilename)
    if not os.path.exists(fastree_treefilename) or os.path.getsize(fastree_treefilename) <= 0:
        # Check whether the msa sliced fasta has enough reads to make a good tree
        if total_slice_seq < 0:
            total_slice_seq = Utility.get_total_seq_from_fasta(msa_window_fasta_filename)
        if total_slice_seq < window_depth_thresh:
            LOGGER.warn("MSA Window " + msa_window_fasta_filename + " does not satisfy window depth constraints")
        else:
            LOGGER.debug("MSA Window " + msa_window_fasta_filename + " satisfies window depth constraints")

            # Feed window fasta into fasttree to make a tree
            os.environ[ENV_OMP_NUM_THREADS] = str(threads)
            with open(fasttree_stdouterr_filename, 'w') as fasttree_stdouterr_fh:
                subprocess.check_call([fastree_exe, '-gtr', '-nt', '-gamma', '-nosupport',
                                       '-log', fastree_logfilename, '-out', fastree_treefilename,
                                       msa_window_fasta_filename],
                                      stdout=fasttree_stdouterr_fh, stderr=fasttree_stdouterr_fh, shell=False,
                                      env=os.environ)

        LOGGER.debug("Done Fasttree for window " + fastree_treefilename)
    else:
        LOGGER.debug("Found existing Fasttree for window " + fastree_treefilename + ". Not regenerating")

    if total_slice_seq >= window_depth_thresh:
        do_hyphy(hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, threads=threads, hyphy_filename_prefix=msa_window_filename_prefix,
             mode=mode, codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename, pvalue=pvalue)


def eval_windows_async(ref, ref_len, sam_filename, out_dir,
                       mapping_cutoff, read_qual_cutoff, max_prop_N,
                       start_nucpos, end_nucpos,
                       windowsize, window_depth_thresh, window_breadth_thresh,
                       pvalue, threads_per_window, concurrent_windows, output_dnds_tsv_filename=None,
                       mode="DNDS", window_slide=3,
                       hyphy_exe=HYPHY_EXE, hyphy_basedir=HYPHY_BASEDIR, fastree_exe=FASTTREE_EXE):
    """
    Launch a separate process to analyze each window.
    Each window can use up to <threads_per_window> threads.

    The maximum number of CPU cores required on a single node = threads_per_window x concurrent_windows.

    :param str ref:  reference name
    :param int ref_len:  length of reference in nucleotide bases
    :param str sam_filename:  full filepath to sam file of read alignments against the reference
    :param str out_dir:  output directory
    :param int mapping_cutoff:  mapping quality threshold below which alignments are thrown out
    :param int read_qual_cutoff:  read quality threshold below which bases are converted to N's
    :param float max_prop_N:  maximum fraction of merged paired-end read that are N's.  Below this threshold the read pair is thrown out.
    :param int start_nucpos:  1-based start nucleotide position of the window
    :param int end_nucpos:  1-based end nucleotide position of the window
    :param int windowsize:  size of each window in nucleotide bases
    :param int window_slide:  Number of nucleotide bases to slide the window by.  Default=3.
    :param int window_depth_thresh:  the minimum number of required reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_thresh:  the minimum fraction of a window that merged paired-end read must cover to be included in the window.
    :param float pvalue:  pvalue threshold for detecting dN/dS (used by HyPhy)
    :param int threads_per_window:  number of threads allotted to processing a single window  (only FastTree and HyPhy will be multithreaded)
    :param int concurrent_windows:  the number of windows to process at the same time.
    :param str output_dnds_tsv_filename:  name of output dN/dS tab separated file generated by HyPhy.  Will be created under out_dir.
    :param str hyphy_exe:  full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe:  full filepath to FastTreeMP executable
    """

    pool = pool_traceback.LoggingPool(processes=concurrent_windows)

    # Create a pseudo multiple-sequence aligned fasta file
    msa_fasta_filename = create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref, ref_len=ref_len,
                                               mapping_cutoff=mapping_cutoff, read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N)

    # All nucleotide positions are 1-based
    last_window_start_nucpos = min(end_nucpos, ref_len - windowsize + 1)
    total_windows = (last_window_start_nucpos - start_nucpos + 1)/Utility.NUC_PER_CODON
    process_results = []
    for start_window_nucpos in range(start_nucpos, last_window_start_nucpos+1, window_slide):
        end_window_nucpos = start_window_nucpos + windowsize - 1
        process_args = (msa_fasta_filename, window_depth_thresh, window_breadth_thresh,
                        start_window_nucpos, end_window_nucpos, pvalue, threads_per_window,
                        mode, hyphy_exe, hyphy_basedir, fastree_exe)
        process_result = pool.apply_async(func=eval_window, args=process_args)
        process_results.append(process_result)

    pool.close()
    LOGGER.debug("Done launching window queue.  Wait for them to finish.")
    pool.join()

    for process_result in process_results:
        try:
            process_result.get()
        except Exception, e:
            LOGGER.error("Error in one of child processes:\n" + e.message)

    LOGGER.debug("Done waiting for window queue.  About to tabulate dn/ds results.")

    if mode == "DNDS":
        dnds_tsv_comments = ("ref=" + ref + ","
                             "ref_len=" + str(ref_len) + "," +
                             "sam=" + sam_filename + "," +
                             "mapping qual cutoff=" + str(mapping_cutoff) + "," +
                             "read qual cutoff=" + str(read_qual_cutoff) + "," +
                             "max fraction N=" + str(max_prop_N) + "," +
                             "start nuc pos=" + str(start_nucpos) + "," +
                             "end nuc pos=" + str(end_nucpos) + "," +
                             "windowsize=" + str(windowsize) + "," +
                             "window depth thresh=" + str(window_depth_thresh) + "," +
                             "window breadth fraction=" + str(window_breadth_thresh) + "," +
                             "pvalue=" + str(pvalue))
        LOGGER.debug("Start Ave Dn/DS for all windows for ref " + ref + " " + output_dnds_tsv_filename)
        seq_dnds_info = slice_miseq.tabulate_dnds(dnds_tsv_dir=out_dir, pvalue_thresh=pvalue, ref=ref, ref_nuc_len=ref_len,
                                      comments=dnds_tsv_comments, output_dnds_tsv_filename=output_dnds_tsv_filename)
        LOGGER.debug("Done Ave Dn/DS for all windows  for ref " + ref + ".  Wrote to " + output_dnds_tsv_filename)
        return seq_dnds_info
    else:
        LOGGER.debug("Done all windows  for ref " + ref)




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", help="full filepath to sam alignment file")
    parser.add_argument("--ref", help="name of reference contig")
    parser.add_argument("--ref_len", type=int, help="length of reference contig in nucleotide bases")
    parser.add_argument("--out_dir", help="output directory in which the pipeline will write all its intermediate files")
    parser.add_argument("--map_q_thresh", type=int, help="mapping quality threshold below which alignments are ignored")
    parser.add_argument("--read_q_thresh", type=int, help="read quality threshold below which bases are converted to Ns")
    parser.add_argument("--max_N_thresh", type=float,
                        help="maximum fraction of Ns allowed in the merged paired-end read below which the paired-end read is ignored")
    parser.add_argument("--window_size", type=int, help="window size in nucleotides")
    parser.add_argument("--window_slide", type=int, help="Number of bases to slide each window by")
    parser.add_argument("--window_breadth_thresh", type=float,
                        help="fraction of window that merged paired-end read must cover with non-gap and non-N nucleotides.  Below this threshold, the read is omitted from the window.")
    parser.add_argument("--window_depth_thresh", type=int,
                        help="1-based start nucleotide position in the reference contig.  The first window will start at this position.")
    parser.add_argument("--start_nucpos", type=int,
                        help="1-based start nucleotide position in the reference contig.  The first window will start at this position.")
    parser.add_argument("--end_nucpos", type=int,
                        help="1-based end nucleotide position in the reference contig.  The last window will start at or before this position.")
    parser.add_argument("--pvalue", type=float,
                        help=" p-value threshold for determining selection significantly different from neutral evolution.")
    parser.add_argument("--threads_per_window", type=int,
                        help="threads allotted per window.  Default: 1")
    parser.add_argument("--concurrent_windows", type=int, help="max number of windows to process concurrently.  Default: 1")
    parser.add_argument("--dnds_tsv",
                        help="full filepath of final tab-separated values file containing selection information for each codon site in the reference from averaged over multiple windows")
    parser.add_argument("--hyphy_exe", help="full filepath of HYPHYMP executable.  Default: taken PATH")
    parser.add_argument("--hyphy_basedir", help="full filepath of HyPhy base directory containing template batch files.  Default: /usr/local/lib/hyphy/TemplateBatchFiles/")
    parser.add_argument("--fastree_exe", help="full filepath of FastTreeMP executable.  Default: taken from PATH")
    parser.add_argument("--mode", help="DNDS or NUC_SUBST.  Default: DNDS")


    args = parser.parse_args()

    eval_windows_args = {}
    if args.ref:
        eval_windows_args["ref"] = args.ref
    if args.ref_len:
        eval_windows_args["ref_len"] = args.ref_len
    if args.sam:
        eval_windows_args["sam_filename"] = args.sam
    if args.out_dir:
        eval_windows_args["out_dir"] = args.out_dir
    if args.map_q_thresh:
        eval_windows_args["mapping_cutoff"] = args.map_q_thresh
    if args.read_q_thresh:
        eval_windows_args["read_qual_cutoff"] = args.read_q_thresh
    if args.max_N_thresh:
        eval_windows_args["max_prop_N"] = args.max_N_thresh
    if args.start_nucpos:
        eval_windows_args["start_nucpos"] = args.start_nucpos
    if args.end_nucpos:
        eval_windows_args["end_nucpos"] = args.end_nucpos
    if args.window_size:
        eval_windows_args["windowsize"] = args.window_size
    if args.window_depth_thresh:
        eval_windows_args["window_depth_thresh"] = args.window_depth_thresh
    if args.window_breadth_thresh:
        eval_windows_args["window_breadth_thresh"] = args.window_breadth_thresh
    if args.pvalue:
        eval_windows_args["pvalue"] = args.pvalue
    if args.threads_per_window:
        eval_windows_args["threads_per_window"] = args.threads_per_window
    if args.concurrent_windows:
        eval_windows_args["concurrent_windows"] = args.concurrent_windows
    if args.dnds_tsv:
        eval_windows_args["output_dnds_tsv_filename"] = args.dnds_tsv
    if args.mode:
        eval_windows_args["mode"] = args.mode
    if args.window_slide:
        eval_windows_args["window_slide"] = args.window_slide
    if args.hyphy_exe:
        eval_windows_args["hyphy_exe"] = args.hyphy_exe
    if args.hyphy_basedir:
        eval_windows_args["hyphy_basedir"] = args.hyphy_basedir
    if args.fastree_exe:
        eval_windows_args["fastree_exe"] = args.fastree_exe


    ordered_args = []
    seq_dnds_info = eval_windows_async(*ordered_args, **eval_windows_args)

if __name__ == '__main__':
    main()
