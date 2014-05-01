import slice_miseq
import sam_handler
import os, sys
import subprocess
import Utility
import logging
import argparse
import pool_traceback

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)

SELECTION_BF = "QuickSelectionDetection.bf"
HYPHY_EXE = "HYPHYMP"
HYPHY_BASEDIR = "/usr/local/lib/hyphy/TemplateBatchFiles/"

FASTTREE_EXE = "FastTreeMP"
ENV_OMP_NUM_THREADS = 'OMP_NUM_THREADS'

NUC_PER_CODON = 3

CMDLINE_OPTIONS = [
    "sam=",
    "ref=",
    "ref_len=",
    "out_dir=",
    "map_q_thresh=",
    "read_q_thresh=",
    "max_N_thresh=",
    "window_size=",
    "window_breadth_thresh=",
    "window_depth_thresh=",
    "start_nucpos=",
    "end_nucpos=",
    "pvalue=",
    "threads_per_window=",
    "concurrent_windows=",
    "dnds_tsv=",
    "hyphy_exe",
    "hyphy_basedir",
    "fastree_exe"]



def create_full_msa_fasta(sam_filename, out_dir, ref, ref_len, mapping_cutoff, read_qual_cutoff, max_prop_N):
    # Create a pseudo multiple-sequence aligned fasta file
    # TODO:  handle indels in multiple sequence align file
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

def eval_windows(ref, ref_len, sam_filename, out_dir, output_dnds_tsv_filename,
                 mapping_cutoff, read_qual_cutoff, max_prop_N,
                 start_nucpos, end_nucpos,
                 window_depth_thresh, window_breadth_thresh, windowsize=300,
                 pvalue=0.05, threads=1, ):
    """
    Slides window along genome.
    Creates the multiple sequence aligned fasta files for the window.
    Feeds the multiple-sequence aligned fasta files to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.
    """
    # TODO:  how to get ORFs?
    # TODO:  do not use kalign2 from PATH
    # TODO:  do not hyphy from PATH

    # Create a pseudo multiple-sequence aligned fasta file
    # TODO:  handle indels in multiple sequence align file
    msa_fasta_filename = create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref, ref_len=ref_len,
                          mapping_cutoff=mapping_cutoff, read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N)

    # All nucleotide positions are 1-based
    last_window_start_nucpos = end_nucpos - windowsize
    for start_window_nucpos in range(start_nucpos, last_window_start_nucpos+1, NUC_PER_CODON):
        end_window_nucpos = start_window_nucpos + windowsize - 1

        eval_window(msa_fasta_filename=msa_fasta_filename, out_dir=out_dir,
                    mapping_cutoff=mapping_cutoff, read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N,
                    window_depth_thresh=window_depth_thresh, window_breadth_thresh=window_breadth_thresh,
                    start_nucpos=start_window_nucpos, end_nucpos=end_window_nucpos,
                    ref=ref, ref_len=ref_len,
                    pvalue=pvalue, threads=threads)

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
    seq_dnds_info = tabulate_dnds(dnds_tsv_dir=out_dir, pvalue_thresh=pvalue, ref=ref, ref_nuc_len=ref_len,
                                  comments=dnds_tsv_comments, output_dnds_tsv_filename=output_dnds_tsv_filename)
    LOGGER.debug("Done Ave Dn/DS for all windows  for ref " + ref + ".  Wrote to " + output_dnds_tsv_filename)
    return seq_dnds_info


def eval_window(msa_fasta_filename, window_depth_thresh, window_breadth_thresh, start_nucpos, end_nucpos, pvalue, threads,
                hyphy_exe=HYPHY_EXE, hyphy_basedir=HYPHY_BASEDIR, fastree_exe=FASTTREE_EXE):
    """
    Slides window along genome.
    Creates the multiple sequence aligned fasta files for the window.
    Feeds the multiple-sequence aligned fasta files to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.
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
            subprocess.check_call([fastree_exe, '-gtr', '-nt', '-nosupport',
                                   '-log', fastree_logfilename, '-out', fastree_treefilename,
                                   msa_window_fasta_filename],
                                  stdout=fasttree_stdouterr_fh, stderr=fasttree_stdouterr_fh, shell=False,
                                  env=os.environ)
        LOGGER.debug("Done Fasttree for window " + fastree_treefilename)
    else:
        LOGGER.debug("Found existing Fasttree for window " + fastree_treefilename + ". Not regenerating")

    hyphy_modelfit_filename = msa_window_filename_prefix + ".nucmodelfit"
    hyphy_dnds_tsv_filename = msa_window_filename_prefix + ".dnds.tsv"
    LOGGER.debug("Start HyPhy for window " + hyphy_dnds_tsv_filename)
    if not os.path.exists(hyphy_dnds_tsv_filename) or os.path.getsize(hyphy_dnds_tsv_filename) <= 0:
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
                                    str(pvalue), # pvalue threshold
                                    "2",    # Export to file
                                    os.path.abspath(hyphy_dnds_tsv_filename),   # dN/dS tsv output file
                                    "2\n"]) # Count approximate numbers of dN, dS rate classes supported by data

        # Feed window tree into hyphy to find dnds for the window
        hyphy_log = msa_window_filename_prefix + ".hyphy.log"
        with open(hyphy_log, 'w') as hyphy_log_fh:
            hyphy_cmd = [hyphy_exe, "BASEPATH="+hyphy_basedir, "CPU=" + str(threads), SELECTION_BF]
            hyphy_proc = subprocess.Popen(hyphy_cmd,
                                          stdin=subprocess.PIPE, stdout=hyphy_log_fh, stderr=hyphy_log_fh, shell=False)
            hyphy_proc.communicate(hyphy_input_str)

            if hyphy_proc.returncode:
                raise subprocess.CalledProcessError(cmd=hyphy_cmd, returncode=hyphy_proc.returncode)
        LOGGER.debug("Done HyPhy for window " + hyphy_dnds_tsv_filename)
    else:
        LOGGER.debug("Found existing HyPhy for window " + hyphy_dnds_tsv_filename + ". Not regenerating")


def eval_windows_async(ref, ref_len, sam_filename, out_dir,
                       mapping_cutoff, read_qual_cutoff, max_prop_N,
                       start_nucpos, end_nucpos,
                       windowsize, window_depth_thresh, window_breadth_thresh,
                       pvalue, threads_per_window, concurrent_windows, output_dnds_tsv_filename,
                       hyphy_exe, hyphy_basedir, fastree_exe):
    """
    Launch a separate process to analyze each window.
    Each window can use up to <threads_per_window> threads.
    """

    pool = pool_traceback.LoggingPool(processes=concurrent_windows)

    # Create a pseudo multiple-sequence aligned fasta file
    # TODO:  handle indels in multiple sequence align file
    msa_fasta_filename = create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref, ref_len=ref_len,
                                                mapping_cutoff=mapping_cutoff, read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N)

    # All nucleotide positions are 1-based
    last_window_start_nucpos = min(end_nucpos, ref_len - windowsize + 1)
    total_windows = (last_window_start_nucpos - start_nucpos + 1)/NUC_PER_CODON
    process_results = []
    for start_window_nucpos in range(start_nucpos, last_window_start_nucpos+1, NUC_PER_CODON):
        end_window_nucpos = start_window_nucpos + windowsize - 1
        process_args = (msa_fasta_filename, window_depth_thresh, window_breadth_thresh,
                        start_window_nucpos, end_window_nucpos, pvalue, threads_per_window,
                        hyphy_exe, hyphy_basedir, fastree_exe)
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam")
    parser.add_argument("--ref")
    parser.add_argument("--ref_len", type=int)
    parser.add_argument("--out_dir")
    parser.add_argument("--map_q_thresh", type=int)
    parser.add_argument("--read_q_thresh", type=int)
    parser.add_argument("--max_N_thresh", type=float)
    parser.add_argument("--window_size", type=int)
    parser.add_argument("--window_breadth_thresh", type=float)
    parser.add_argument("--window_depth_thresh", type=int)
    parser.add_argument("--start_nucpos", type=int)
    parser.add_argument("--end_nucpos", type=int)
    parser.add_argument("--pvalue", type=float)
    parser.add_argument("--threads_per_window", type=int)
    parser.add_argument("--concurrent_windows", type=int)
    parser.add_argument("--dnds_tsv")
    parser.add_argument("--hyphy_exe")
    parser.add_argument("--hyphy_basedir")
    parser.add_argument("--fastree_exe")

    args = parser.parse_args()
    print args

    seq_dnds_info = eval_windows_async(sam_filename=args.sam,
                                       ref=args.ref,
                                       ref_len=args.ref_len,
                                       out_dir=args.out_dir,
                                       mapping_cutoff=args.map_q_thresh,
                                       read_qual_cutoff=args.read_q_thresh,
                                       max_prop_N=args.max_N_thresh,
                                       windowsize=args.window_size,
                                       window_breadth_thresh=args.window_breadth_thresh,
                                       window_depth_thresh=args.window_depth_thresh,
                                       start_nucpos=args.start_nucpos,
                                       end_nucpos=args.end_nucpos,
                                       pvalue=args.pvalue,
                                       threads_per_window=args.threads_per_window,
                                       concurrent_windows=args.concurrent_windows,
                                       output_dnds_tsv_filename=args.dnds_tsv,
                                       hyphy_exe=args.hyphy_exe,
                                       hyphy_basedir=args.hyphy_basedir,
                                       fastree_exe=args.fastree_exe)

if __name__ == '__main__':
    main()
