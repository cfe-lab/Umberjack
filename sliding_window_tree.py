import slice_miseq
import sam_handler
import os, sys
import subprocess
import Utility
import logging
import argparse
import pool_traceback
import traceback
import hyphy.hyphy_handler as hyphy
import fasttree.fasttree_handler as fasttree


from mpi4py import MPI


import array


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)


# For MPI
MASTER_RANK = 0
TAG_WORK = 1
TAG_DIE = 2

MODE_DNDS = "DNDS"
MODE_GTR_CMP = "GTR_CMP"
MODE_GTR_RATE = "GTR_RATE"

# TODO:  handle inserts in multiple sequence align file
# TODO:  handle reads that extend the reference in MSA file
# TODO: this is a performance killing step.  Instead of writing all MSA aligned reads to file, we need to just do it for
# parts of the genome that are used in the sliding window or hold all of this in memory.
def create_full_msa_fasta(sam_filename, out_dir, ref, ref_len, mapping_cutoff, read_qual_cutoff, max_prop_N):
    """
    Creates a pseudo multiple-sequence aligned fasta file for all reads using pairwise alignment from a SAM file.


    :return: path to multiple sequence aligned fasta file of all reads
    :rtype : str
    :param str sam_filename: filepath to sam file of read alignments to the reference
    :param str out_dir: output directory
    :param str ref: reference name
    :param int ref_len: length of reference in nucleotides
    :param int mapping_cutoff: mapping quality cutoff
    :param int read_qual_cutoff: read quality cutoff.  Bases below this cutoff are converted to N's.
    :param float max_prop_N:  maximum fraction of bases in a merged mate-pair read that is allowed to be N's.
                                Reads that exceed this threshold are thrown out.
    """

    sam_filename_nopath = os.path.split(sam_filename)[1]
    sam_filename_prefix = os.path.splitext(sam_filename_nopath)[0]
    msa_fasta_filename = out_dir + os.sep + sam_filename_prefix + "." + ref + ".msa.fasta"

    LOGGER.debug("Start Full MSA-Fasta from SAM for ref " + ref)
    if not os.path.exists(msa_fasta_filename) or os.path.getsize(msa_fasta_filename) <= 0:
        sam_handler.create_msa_fasta_from_sam(sam_filename=sam_filename, ref=ref, ref_len=ref_len,
                                              out_fasta_filename=msa_fasta_filename, mapping_cutoff=mapping_cutoff,
                                              read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N)
        LOGGER.debug("Done Full MSA-Fasta from SAM for ref " + ref)
    else:
        LOGGER.warn("Found existing Full MSA-Fasta from SAM for ref " + ref + ".  Not regenerating")

    return msa_fasta_filename



# TODO:  do multiple test corrections for pvalues
def eval_window(msa_fasta_filename, window_depth_cutoff, window_breadth_cutoff, start_window_nucpos, end_window_nucpos,
                pvalue, threads_per_window, mode="DNDS", hyphy_exe=hyphy.HYPHY_EXE, hyphy_basedir=hyphy.HYPHY_BASEDIR,
                fastree_exe=fasttree.FASTTREE_EXE):
    """
    Handles the processing for a single window along the genome.
    Creates the multiple sequence aligned fasta file for the window.
    Feeds the window multiple-sequence aligned fasta file to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.

    :param str msa_fasta_filename: full filepath to multiple sequence aligned file for all reads.
    :param int window_depth_cutoff:  the minimum number of required reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_cutoff: the minimum fraction of a window that merged paired-end read must cover to be included in the window.
    :param int start_window_nucpos:  1-based start nucleotide position of the window
    :param int end_window_nucpos:  1-based end nucleotide position of the window
    :param float pvalue:  pvalue threshold for detecting dN/dS (used by HyPhy)
    :param int threads_per_window: number of threads allotted to processing this window  (only FastTree and HyPhy will be multithreaded)
    :param str hyphy_exe: full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe: full filepath to FastTreeMP executable
    """

    LOGGER.debug("msa_fasta_filename=" + msa_fasta_filename + "\n" +
                 "window_depth_thresh=" + str(window_depth_cutoff) + "\n" +
                 "window_breadth_thresh=" + str(window_breadth_cutoff) + "\n" +
                 "start_nucpos=" + str(start_window_nucpos) + "\n" +
                 "end_nucpos=" + str(end_window_nucpos) + "\n" +
                 "pvalue=" + str(pvalue) + "\n" +
                 "threads=" + str(threads_per_window) + "\n")

    # Slice the multiple sequence aligned fasta file into a window fasta
    msa_fasta_filename_prefix = os.path.splitext(msa_fasta_filename)[0]
    msa_window_filename_prefix = msa_fasta_filename_prefix + "." + str(start_window_nucpos) + "_" + str(end_window_nucpos)
    msa_window_fasta_filename = msa_window_filename_prefix + ".fasta"
    total_slice_seq = slice_miseq.create_slice_msa_fasta(fasta_filename=msa_fasta_filename,
                                                             out_fasta_filename=msa_window_fasta_filename,
                                                             start_pos=start_window_nucpos, end_pos=end_window_nucpos,
                                                             breadth_thresh=window_breadth_cutoff)

    # Check whether the msa sliced fasta has enough reads to make a good tree
    if total_slice_seq < window_depth_cutoff:
        LOGGER.warn("MSA Window " + msa_window_fasta_filename + " does not satisfy window depth constraints")
        return

    LOGGER.debug("MSA Window " + msa_window_fasta_filename + " satisfies window depth constraints")

    fastree_treefilename = fasttree.make_tree(fasta_fname=msa_window_fasta_filename, threads=threads_per_window, fastree_exe=fastree_exe)

    if mode == MODE_DNDS:
        hyphy.calc_dnds(hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, threads=threads_per_window,
                        codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename, pvalue=pvalue)
    elif mode == MODE_GTR_CMP:
        hyphy.calc_nuc_subst(hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, threads=threads_per_window,
                             codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename)
    elif mode != MODE_GTR_RATE:
        raise  ValueError("Invalid mode " + mode)



# TODO:  clean my parameters
def tabulate_results(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                     end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, pvalue, output_csv_filename,
                     mode, window_slide, smooth_dist):

    comments = ("ref=" + ref + ","
                             "ref_len=" + str(ref_len) + "," +
                             "sam=" + sam_filename + "," +
                             "map_qual_cutoff=" + str(map_qual_cutoff) + "," +
                             "read qual cutoff=" + str(read_qual_cutoff) + "," +
                             "max_prop_n=" + str(max_prop_n) + "," +
                             "start nuc pos=" + str(start_nucpos) + "," +
                             "end nuc pos=" + str(end_nucpos) + "," +
                             "windowsize=" + str(window_size) + "," +
                             "window_slide=" + str(window_slide) + "," +
                             "window_depth_cutoff=" + str(window_depth_cutoff) + "," +
                             "window_breadth_cutoff=" + str(window_breadth_cutoff) + "," +
                             "pvalue=" + str(pvalue))
    if mode == MODE_DNDS:
        LOGGER.debug("Start Ave Dn/DS for all windows for ref " + ref + " " + output_csv_filename)
        seq_dnds_info = slice_miseq.tabulate_dnds(dnds_tsv_dir=out_dir, ref=ref, ref_nuc_len=ref_len,
                                                  pvalue_thresh=pvalue, output_csv_filename=output_csv_filename,
                                                  comments=comments, smooth_dist=smooth_dist)
        LOGGER.debug("Done Ave Dn/DS for all windows  for ref " + ref + ".  Wrote to " + output_csv_filename)
        return seq_dnds_info
    elif mode == MODE_GTR_RATE:
        LOGGER.debug("Start Tabulate GTR Rates for All Windows For Ref " + ref + " " + output_csv_filename)
        slice_miseq.tabulate_rates(fasttree_output_dir=out_dir, output_csv_filename=output_csv_filename, comments=comments)
        LOGGER.debug("Done Tabulate GTR Rates for all windows for ref " + ref + " " + output_csv_filename)
    else:
        LOGGER.debug("Done all windows  for ref " + ref)


def eval_windows_async(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                       end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, pvalue, threads_per_window,
                       concurrent_windows, output_csv_filename=None, mode="DNDS", window_slide=3, smooth_dist=50,
                       hyphy_exe=hyphy.HYPHY_EXE, hyphy_basedir=hyphy.HYPHY_BASEDIR, fastree_exe=fasttree.FASTTREE_EXE):
    """
    Launch a separate process to analyze each window.
    Each window can use up to <threads_per_window> threads.

    The maximum number of CPU cores required on a single node = threads_per_window x concurrent_windows.

    :param str ref:  reference name
    :param int ref_len:  length of reference in nucleotide bases
    :param str sam_filename:  full filepath to sam file of read alignments against the reference
    :param str out_dir:  output directory
    :param int map_qual_cutoff:  mapping quality threshold below which alignments are thrown out
    :param int read_qual_cutoff:  read quality threshold below which bases are converted to N's
    :param float max_prop_n:  maximum fraction of merged paired-end read that are N's.  Below this threshold the read pair is thrown out.
    :param int start_nucpos:  1-based start nucleotide position of the window
    :param int end_nucpos:  1-based end nucleotide position of the window
    :param int window_size:  size of each window in nucleotide bases
    :param int window_slide:  Number of nucleotide bases to slide the window by.  Default=3.
    :param int window_depth_cutoff:  the minimum number of required reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_cutoff:  the minimum fraction of a window that merged paired-end read must cover to be included in the window.
    :param float pvalue:  pvalue threshold for detecting dN/dS (used by HyPhy)
    :param int threads_per_window:  number of threads allotted to processing a single window  (only FastTree and HyPhy will be multithreaded)
    :param int concurrent_windows:  the number of windows to process at the same time.
    :param str output_csv_filename:  name of output dN/dS tab separated file generated by HyPhy.  Will be created under out_dir.
    :param str hyphy_exe:  full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe:  full filepath to FastTreeMP executable
    """



    pool = pool_traceback.LoggingPool(processes=concurrent_windows)

    # Create a pseudo multiple-sequence aligned fasta file
    msa_fasta_filename = create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref, ref_len=ref_len,
                                               mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
                                               max_prop_N=max_prop_n)

    # All nucleotide positions are 1-based
    last_window_start_nucpos = min(end_nucpos, ref_len - window_size + 1)

    total_windows = (last_window_start_nucpos - start_nucpos + 1) / Utility.NUC_PER_CODON
    LOGGER.debug("There are " + str(total_windows) + " total windows to process")

    process_results = []
    for start_window_nucpos in range(start_nucpos, last_window_start_nucpos + 1, window_slide):
        end_window_nucpos = start_window_nucpos + window_size - 1
        window_args = {"msa_fasta_filename": msa_fasta_filename,
                       "window_depth_cutoff": window_depth_cutoff,
                       "window_breadth_cutoff": window_breadth_cutoff,
                       "start_window_nucpos": start_window_nucpos,
                       "end_window_nucpos": end_window_nucpos,
                       "pvalue": pvalue,
                       "threads_per_window": threads_per_window,
                       "mode": mode,
                       "hyphy_exe": hyphy_exe,
                       "hyphy_basedir": hyphy_basedir,
                       "fastree_exe": fastree_exe}
        process_result = pool.apply_async(eval_window, (), window_args)
        process_results.append(process_result)

    pool.close()
    LOGGER.debug("Done launching window queue.  Wait for them to finish.")
    pool.join()

    for process_result in process_results:
        try:
            process_result.get()
        except Exception, e:
            LOGGER.error("Error in one of slave processes:\n" + e.message)

    LOGGER.debug("Done waiting for window queue.  About to tabulate results.")

    tabulate_results(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                     end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, pvalue,
                     output_csv_filename, mode, window_slide, smooth_dist)


class WindowSlaveInfo:
    """
    Keeps track of the slave information
    """
    def __init__(self, slave_rank, work_arguments, mpi_send_request, mpi_rcv_request):
        """
        :param slave_rank: integer slave rank (starts from 1)
        :param dict work_arguments:  dict of arguments sent to the slave to do work
        :param mpi_send_request: mpi send request sent to the slave
        :param mpi_rcv_request: mpi receive request sent from the slave
        """
        self.slave_rank = slave_rank
        self.work_arguments = work_arguments
        self.mpi_send_request = mpi_send_request
        self.mpi_rcv_request = mpi_rcv_request

def eval_windows_mpi(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                     end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, pvalue, threads_per_window,
                     output_csv_filename=None, mode="DNDS", window_slide=3, smooth_dist=50, hyphy_exe=hyphy.HYPHY_EXE,
                     hyphy_basedir=hyphy.HYPHY_BASEDIR, fastree_exe=fasttree.FASTTREE_EXE):
    """
    Launch a separate process to analyze each window via MPI.  Similar to eval_windows_async, but uses MPI.

    :param str ref:  reference name
    :param int ref_len:  length of reference in nucleotide bases
    :param str sam_filename:  full filepath to sam file of read alignments against the reference
    :param str out_dir:  output directory
    :param int map_qual_cutoff:  mapping quality threshold below which alignments are thrown out
    :param int read_qual_cutoff:  read quality threshold below which bases are converted to N's
    :param float max_prop_n:  maximum fraction of merged paired-end read that are N's.  Below this threshold the read pair is thrown out.
    :param int start_nucpos:  1-based start nucleotide position of the window
    :param int end_nucpos:  1-based end nucleotide position of the window
    :param int window_size:  size of each window in nucleotide bases
    :param int window_slide:  Number of nucleotide bases to slide the window by.  Default=3.
    :param int window_depth_cutoff:  the minimum number of required reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_cutoff:  the minimum fraction of a window that merged paired-end read must cover to be included in the window.
    :param float pvalue:  pvalue threshold for detecting dN/dS (used by HyPhy)
    :param int threads_per_window:  number of threads allotted to processing a single window  (only FastTree and HyPhy will be multithreaded)
    :param str output_csv_filename:  name of output dN/dS tab separated file generated by HyPhy.  Will be created under out_dir.
    :param str hyphy_exe:  full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe:  full filepath to FastTreeMP executable
    """

    # from mpi4py import MPI

    try:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        LOGGER.debug("I am rank=" + str(rank))
        LOGGER.debug("I am on machine=" + str(MPI.Get_processor_name()))

        if rank == MASTER_RANK:
            pool_size = comm.Get_size()
            LOGGER.debug("Pool size = " + str(pool_size))

            # Create a pseudo multiple-sequence aligned fasta file
            msa_fasta_filename = create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref, ref_len=ref_len,
                                                       mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
                                                       max_prop_N=max_prop_n)

            # All nucleotide positions are 1-based
            last_window_start_nucpos = min(end_nucpos, ref_len - window_size + 1)

            total_windows = (last_window_start_nucpos - start_nucpos + 1) / Utility.NUC_PER_CODON
            LOGGER.debug("Launching " + str(total_windows) + " total windows")

            available_slaves = range(1, pool_size)
            busy_slave_2_request = {}

            start_window_nucpos = start_nucpos



            while start_window_nucpos <= last_window_start_nucpos or busy_slave_2_request:

                # Assign work to slaves
                while start_window_nucpos <= last_window_start_nucpos and available_slaves:
                    end_window_nucpos = start_window_nucpos + window_size - 1

                    window_args = {"msa_fasta_filename": msa_fasta_filename,
                                   "window_depth_cutoff": window_depth_cutoff,
                                   "window_breadth_cutoff": window_breadth_cutoff,
                                   "start_window_nucpos": start_window_nucpos,
                                   "end_window_nucpos": end_window_nucpos,
                                   "pvalue": pvalue,
                                   "threads_per_window": threads_per_window,
                                   "mode": mode,
                                   "hyphy_exe": hyphy_exe,
                                   "hyphy_basedir": hyphy_basedir,
                                   "fastree_exe": fastree_exe}
                    slave_rank = available_slaves.pop(0)

                    str_window_args = ', '.join('{}:{}'.format(key, val) for key, val in window_args.items())
                    LOGGER.debug("Sending window_args to slave=" + str(slave_rank) + " " + str_window_args)

                    send_request = comm.isend(obj=window_args, dest=slave_rank, tag=TAG_WORK)  # non-blocking
                    rcv_request = comm.irecv(dest=slave_rank, tag=MPI.ANY_TAG)  # non-blocking

                    # MPI.Comm.isend() stores a copy of the pickled (i.e. serialized) windows_args dict
                    # in the buffer of a new MPI.Request object  (send_request)
                    # The slave accesses the buffer in the send_request MPI.Request object when it calls MPI.Comm.recv()
                    # to retrieve work.
                    # If we do not keep send_request in scope, its buffer can be overwitten by a random process by the time
                    # the slave gets to it.
                    # When the master does a non-blocking call to retrieve an response form the slave in MPI.Comm.irecv,
                    # it gets another MPI.Request object (rcv_request) whose buffer will contain the response from the slave when it's finished.
                    # The buffers from send_request and rcv_request are completely separate from each other.
                    busy_slave_2_request[slave_rank] = WindowSlaveInfo(slave_rank=slave_rank,
                                                                       work_arguments=window_args,
                                                                       mpi_send_request=send_request,
                                                                       mpi_rcv_request=rcv_request)

                    start_window_nucpos += window_slide

                # Check on slaves
                if busy_slave_2_request:
                    requests = [window_slave_info.mpi_rcv_request for window_slave_info in busy_slave_2_request.values()]
                    mpi_status = MPI.Status()
                    idx, err_msg = MPI.Request.waitany(requests=requests, status=mpi_status)
                    done_slave_rank = mpi_status.Get_source()
                    available_slaves.extend([done_slave_rank])
                    busy_slave_2_request.pop(done_slave_rank)
                    if err_msg:
                        LOGGER.error("Received error from slave=" + str(done_slave_rank) + " err_msg=" + str(err_msg))
                    else:
                        LOGGER.debug("Received success from slave=" + str(done_slave_rank))

            LOGGER.debug("Done Launching " + str(total_windows) + " total windows")

            LOGGER.debug("About to Kill Slaves")
            for slave_rank in range(1, pool_size):
                comm.isend(obj=None, dest=slave_rank, tag=TAG_DIE)
            LOGGER.debug("Done Killing slaves.")

            LOGGER.debug("About to tabulate results")
            tabulate_results(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n,
                             start_nucpos, end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, pvalue,
                             output_csv_filename, mode, window_slide, smooth_dist)
            LOGGER.debug("Done tabulating results")

        else:  # slave process does the work
            is_die = False
            while not is_die:
                try:
                    mpi_status = MPI.Status()
                    window_args = comm.recv(source=MASTER_RANK, tag=MPI.ANY_TAG, status=mpi_status)  # block till the master tells me to do something


                    if mpi_status.Get_tag() == TAG_DIE:  # master wants me to die
                        is_die = True
                        LOGGER.debug("Master wants me to die - I was rank " + str(rank))
                    else:  # master wants me to work
                        str_window_args = ', '.join('{}:{}'.format(key, val) for key, val in window_args.items())
                        LOGGER.debug("Received window_args=" + str_window_args)
                        eval_window(**window_args)
                        comm.send(obj=None, dest=MASTER_RANK, tag=TAG_WORK)  # Tell master that I'm done and wait for master to receive my response
                except Exception, e:
                    LOGGER.exception("Failure in slave=" + str(rank))
                    err_msg = traceback.format_exc()
                    comm.send(obj=err_msg, dest=MASTER_RANK, tag=TAG_WORK)  # Tell master that I encountered an exception
    except Exception, e:
        LOGGER.exception("Uncaught Exception.  Aborting")
        comm.Abort()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam_filename", help="full filepath to sam alignment file")
    parser.add_argument("--ref", help="name of reference contig")
    parser.add_argument("--ref_len", type=int, help="length of reference contig in nucleotide bases")
    parser.add_argument("--out_dir",
                        help="output directory in which the pipeline will write all its intermediate files")
    parser.add_argument("--map_qual_cutoff", type=int,
                        help="mapping quality threshold below which alignments are ignored")
    parser.add_argument("--read_qual_cutoff", type=int,
                        help="read quality threshold below which bases are converted to Ns")
    parser.add_argument("--max_prop_n", type=float,
                        help="maximum fraction of Ns allowed in the merged paired-end read below which the paired-end read is ignored")
    parser.add_argument("--window_size", type=int, help="window size in nucleotides")
    parser.add_argument("--window_slide", type=int, help="Number of bases to slide each window by")
    parser.add_argument("--smooth_dist", type=int, help="Site dn/ds is smoothed by this many codons on either side")
    parser.add_argument("--window_breadth_cutoff", type=float,
                        help="fraction of window that merged paired-end read must cover with non-gap and non-N nucleotides.  Below this threshold, the read is omitted from the window.")
    parser.add_argument("--window_depth_cutoff", type=int,
                        help="1-based start nucleotide position in the reference contig.  The first window will start at this position.")
    parser.add_argument("--start_nucpos", type=int,
                        help="1-based start nucleotide position in the reference contig.  The first window will start at this position.")
    parser.add_argument("--end_nucpos", type=int,
                        help="1-based end nucleotide position in the reference contig.  The last window will start at or before this position.")
    parser.add_argument("--pvalue", type=float,
                        help=" p-value threshold for determining selection significantly different from neutral evolution.")
    parser.add_argument("--threads_per_window", type=int,
                        help="threads allotted per window.  Default: 1")
    parser.add_argument("--concurrent_windows", type=int, help="Max number of windows to process concurrently.  " +
                                                               "Ignored when --mpi is defined.  Default: 1")
    parser.add_argument("--output_csv_filename",
                        help="full filepath of final comma-separated values file containing selection information for each codon site in the reference from averaged over multiple windows")
    parser.add_argument("--hyphy_exe", help="full filepath of HYPHYMP executable.  Default: taken PATH")
    parser.add_argument("--hyphy_basedir",
                        help="full filepath of HyPhy base directory containing template batch files.  Default: /usr/local/lib/hyphy/TemplateBatchFiles/")
    parser.add_argument("--fastree_exe", help="full filepath of FastTreeMP executable.  Default: taken from PATH")
    parser.add_argument("--mode", help="DNDS or GTR_CMP.  Default: DNDS")
    parser.add_argument("--mpi", help="Runs in MPI mode with multiple processes on multiple nodes." +
                                      "If python module mpi4py is not installed, then runs multiple processes on single node. " +
                                      "Default: false")


    args = parser.parse_args()
    eval_windows_args = {}
    for key, val in vars(args).iteritems():
        if val:
            eval_windows_args[key] = val

    # if the user has mpi4py installed, then runs the MPI version
    # otherwise runs the multiprocessing version on current node
    do_mpi = False
    if args.mpi:
        try:
            from mpi4py import MPI
            # Ignore the concurrent_windows commandline arg and uses the number of processors indicated by mpirun command
            eval_windows_args.pop("concurrent_windows", None)
            LOGGER.debug("Running MPI Version. Ignoring --concurrent_windows flag.  Using mpirun node arguments.")
            do_mpi = True
        except ImportError:
            LOGGER.warn("You must install mpi4py module in order to leverage multiple nodes.  Running on single node.")

    eval_windows_args.pop("mpi", None)  # this is not used in eval_windows* methods
    if do_mpi:
        eval_windows_mpi(**eval_windows_args)
    else:
        eval_windows_async(**eval_windows_args)


if __name__ == '__main__':
    main()
