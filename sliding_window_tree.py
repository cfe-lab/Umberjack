import os
import sys
import logging
import argparse
import traceback

from sam import sam_handler
import slice_miseq
import Utility
import pool_traceback
import hyphy.hyphy_handler as hyphy
import fasttree.fasttree_handler as fasttree


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)


# For MPI
PRIMARY_RANK = 0
TAG_WORK = 1
TAG_TERMINATE = 2

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
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
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



def eval_window(window_depth_cutoff, window_breadth_cutoff, start_window_nucpos, end_window_nucpos, ref, out_dir, pvalue=0.05,
                msa_fasta_filename=None, sam_filename=None, map_qual_cutoff=None, read_qual_cutoff=None, max_prop_N=None, threads_per_window=1, mode="DNDS",
                hyphy_exe=hyphy.HYPHY_EXE, hyphy_basedir=hyphy.HYPHY_BASEDIR, fastree_exe=fasttree.FASTTREEMP_EXE):
    """
    Handles the processing for a single window along the genome.
    Creates the multiple sequence aligned fasta file for the window.
    Feeds the window multiple-sequence aligned fasta file to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.

    :param sam_filename:
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

    LOGGER.debug("msa_fasta_filename=" + str(msa_fasta_filename) + "\n" +
                 "sam_filename=" + str(sam_filename) + "\n" +
                 "window_depth_thresh=" + str(window_depth_cutoff) + "\n" +
                 "window_breadth_thresh=" + str(window_breadth_cutoff) + "\n" +
                 "start_nucpos=" + str(start_window_nucpos) + "\n" +
                 "end_nucpos=" + str(end_window_nucpos) + "\n" +
                 "pvalue=" + str(pvalue) + "\n" +
                 "threads=" + str(threads_per_window) + "\n")

    if not msa_fasta_filename and not sam_filename:
        raise ValueError("Either msa_fasta_filename or sam_filename must be defined")

    if sam_filename:
        sam_filename_nopath = os.path.split(sam_filename)[1]
        sam_filename_prefix = os.path.splitext(sam_filename_nopath)[0]

        msa_window_filename_prefix = out_dir + os.sep + sam_filename_prefix + "." + str(start_window_nucpos) + "_" + str(end_window_nucpos)
        msa_window_fasta_filename = msa_window_filename_prefix + ".fasta"
        total_slice_seq = sam_handler.create_msa_slice_from_sam(sam_filename=sam_filename,
                                                                ref=ref,
                                                                out_fasta_filename=msa_window_fasta_filename,
                                                                mapping_cutoff=map_qual_cutoff,
                                                                read_qual_cutoff=read_qual_cutoff,
                                                                max_prop_N=max_prop_N,
                                                                breadth_thresh=window_breadth_cutoff,
                                                                start_pos=start_window_nucpos,
                                                                end_pos=end_window_nucpos,
                                                                is_insert=False)
    else:
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
        hyphy.calc_dnds(codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename,
                        hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, threads=threads_per_window)
    elif mode == MODE_GTR_CMP:
        hyphy.calc_nuc_subst(hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, threads=threads_per_window,
                             codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename)
    elif mode != MODE_GTR_RATE:
        raise  ValueError("Invalid mode " + mode)



# TODO:  clean my parameters
def tabulate_results(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                     end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, output_csv_filename, mode,
                     window_slide, smooth_dist):

    comments = ("ref=" + ref + "," +
                "ref_len=" + str(ref_len) + "," +
                "sam=" + sam_filename + "," +
                "map_qual_cutoff=" + str(map_qual_cutoff) + "," +
                "read qual cutoff=" + str(read_qual_cutoff) + "," +
                "max_prop_n=" + str(max_prop_n) + "," +
                "start_nuc_pos=" + str(start_nucpos) + "," +
                "end_nuc_pos=" + str(end_nucpos) + "," +
                "windowsize=" + str(window_size) + "," +
                "window_slide=" + str(window_slide) + "," +
                "window_depth_cutoff=" + str(window_depth_cutoff) + "," +
                "window_breadth_cutoff=" + str(window_breadth_cutoff) + "," +
                "smooth_dist=" + str(smooth_dist))
    if mode == MODE_DNDS:
        LOGGER.debug("Start Ave Dn/DS for all windows for ref " + ref + " " + output_csv_filename)
        seq_dnds_info = slice_miseq.tabulate_dnds(dnds_tsv_dir=out_dir, ref=ref, ref_nuc_len=ref_len,
                                                  output_csv_filename=output_csv_filename, comments=comments,
                                                  smooth_dist=smooth_dist)
        LOGGER.debug("Done Ave Dn/DS for all windows  for ref " + ref + ".  Wrote to " + output_csv_filename)
        return seq_dnds_info
    elif mode == MODE_GTR_RATE:
        LOGGER.debug("Start Tabulate GTR Rates for All Windows For Ref " + ref + " " + output_csv_filename)
        slice_miseq.tabulate_rates(fasttree_output_dir=out_dir, output_csv_filename=output_csv_filename, comments=comments)
        LOGGER.debug("Done Tabulate GTR Rates for all windows for ref " + ref + " " + output_csv_filename)
    else:
        LOGGER.debug("Done all windows  for ref " + ref)


def eval_windows_async(ref, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                       end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, threads_per_window,
                       concurrent_windows, output_csv_filename=None, mode="DNDS", window_slide=3, smooth_dist=10,
                       hyphy_exe=hyphy.HYPHY_EXE, hyphy_basedir=hyphy.HYPHY_BASEDIR, fastree_exe=fasttree.FASTTREEMP_EXE):
    """
    Launch a separate process to analyze each window.
    Each window can use up to <threads_per_window> threads.

    The maximum number of CPU cores required on a single node = threads_per_window x concurrent_windows.

    :param str ref:  reference name
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
    :param int threads_per_window:  number of threads allotted to processing a single window  (only FastTree and HyPhy will be multithreaded)
    :param int concurrent_windows:  the number of windows to process at the same time.
    :param str output_csv_filename:  name of output dN/dS tab separated file generated by HyPhy.  Will be created under out_dir.
    :param str hyphy_exe:  full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe:  full filepath to FastTreeMP executable
    """



    pool = pool_traceback.LoggingPool(processes=concurrent_windows)

    ref_len = sam_handler.get_ref_len(sam_filename, ref)

    # Create a pseudo multiple-sequence aligned fasta file
    # msa_fasta_filename = create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref, ref_len=ref_len,
    #                                            mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
    #                                            max_prop_N=max_prop_n)

    msa_fasta_filename = None  # TODO:  write the full msa if in debug mode
    # All nucleotide positions are 1-based
    last_window_start_nucpos = min(end_nucpos, ref_len - window_size + 1)

    total_windows = (last_window_start_nucpos - start_nucpos + 1) / Utility.NUC_PER_CODON
    LOGGER.debug("There are " + str(total_windows) + " total windows to process")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)  # TODO:  put msa_fasta_file into the outdir too
    process_results = []
    for start_window_nucpos in range(start_nucpos, last_window_start_nucpos + 1, window_slide):
        end_window_nucpos = start_window_nucpos + window_size - 1
        window_args = {"window_depth_cutoff": window_depth_cutoff,
                       "window_breadth_cutoff": window_breadth_cutoff,
                       "start_window_nucpos": start_window_nucpos,
                       "end_window_nucpos": end_window_nucpos,
                       "ref": ref,
                       "out_dir": out_dir,
                       "msa_fasta_filename": msa_fasta_filename,
                       "sam_filename": sam_filename,
                       "map_qual_cutoff": map_qual_cutoff,
                       "read_qual_cutoff": read_qual_cutoff,
                       "max_prop_N": max_prop_n,
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
            LOGGER.error("Error in one of replica processes:\n" + e.message)

    LOGGER.debug("Done waiting for window queue.  About to tabulate results.")

    tabulate_results(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                     end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, output_csv_filename, mode, window_slide,
                     smooth_dist)



class WindowReplicaInfo:
    """
    Keeps track of the replica information
    """
    def __init__(self, replica_rank, work_arguments, mpi_send_request, mpi_rcv_request):
        """
        :param replica_rank: integer replica rank (starts from 1)
        :param dict work_arguments:  dict of arguments sent to the replica to do work
        :param mpi_send_request: mpi request sent to replica
        :param mpi_rcv_request: mpi request received from replica
        """
        self.replica_rank = replica_rank
        self.work_arguments = work_arguments
        self.mpi_send_request = mpi_send_request
        self.mpi_rcv_request = mpi_rcv_request


def eval_windows_mpi(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                     end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, threads_per_window,
                     output_csv_filename=None, mode="DNDS", window_slide=3, smooth_dist=10, hyphy_exe=hyphy.HYPHY_EXE,
                     hyphy_basedir=hyphy.HYPHY_BASEDIR, fastree_exe=fasttree.FASTTREEMP_EXE):
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
    :param int threads_per_window:  number of threads allotted to processing a single window  (only FastTree and HyPhy will be multithreaded)
    :param str output_csv_filename:  name of output dN/dS tab separated file generated by HyPhy.  Will be created under out_dir.
    :param str hyphy_exe:  full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe:  full filepath to FastTreeMP executable
    """

    from mpi4py import MPI  # TODO:  remove this before checking in
    comm = MPI.COMM_WORLD
    try:

        rank = comm.Get_rank()
        LOGGER.debug("I am rank=" + str(rank))
        LOGGER.debug("I am on machine=" + str(MPI.Get_processor_name()))

        if rank == PRIMARY_RANK:
            pool_size = comm.Get_size()
            LOGGER.debug("Pool size = " + str(pool_size))

            # # Create a pseudo multiple-sequence aligned fasta file
            # msa_fasta_filename = create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref, ref_len=ref_len,
            #                                            mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
            #                                            max_prop_N=max_prop_n)
            msa_fasta_filename = None  # TODO:  only write this out if debug mode
            # All nucleotide positions are 1-based
            last_window_start_nucpos = min(end_nucpos, ref_len - window_size + 1)

            total_windows = (last_window_start_nucpos - start_nucpos + 1) / Utility.NUC_PER_CODON
            LOGGER.debug("Launching " + str(total_windows) + " total windows")

            available_replicas = range(1, pool_size)
            busy_replica_2_request = {}

            start_window_nucpos = start_nucpos

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)  # TODO:  put msa_fasta_file into the outdir too

            while start_window_nucpos <= last_window_start_nucpos or busy_replica_2_request:

                # Assign work to replicas
                while start_window_nucpos <= last_window_start_nucpos and available_replicas:
                    end_window_nucpos = start_window_nucpos + window_size - 1

                    window_args = {"window_depth_cutoff": window_depth_cutoff,
                                   "window_breadth_cutoff": window_breadth_cutoff,
                                   "start_window_nucpos": start_window_nucpos,
                                   "end_window_nucpos": end_window_nucpos,
                                   "ref": ref,
                                   "out_dir": out_dir,
                                   "msa_fasta_filename": msa_fasta_filename,
                                   "sam_filename": sam_filename,
                                   "map_qual_cutoff": map_qual_cutoff,
                                   "read_qual_cutoff": read_qual_cutoff,
                                   "max_prop_N": max_prop_n,
                                   "threads_per_window": threads_per_window,
                                   "mode": mode,
                                   "hyphy_exe": hyphy_exe,
                                   "hyphy_basedir": hyphy_basedir,
                                   "fastree_exe": fastree_exe}
                    replica_rank = available_replicas.pop(0)

                    str_window_args = ', '.join('{}:{}'.format(key, val) for key, val in window_args.items())
                    LOGGER.debug("Sending window_args to replica=" + str(replica_rank) + " " + str_window_args)


                    send_request = comm.isend(obj=window_args, dest=replica_rank, tag=TAG_WORK)  # non-blocking
                    rcv_request = comm.irecv(dest=replica_rank, tag=MPI.ANY_TAG)  # non-blocking

                    # MPI.Comm.isend() stores a copy of the pickled (i.e. serialized) windows_args dict
                    # in the buffer of a new MPI.Request object  (send_request)
                    # The replica accesses the buffer in the send_request MPI.Request object when it calls MPI.Comm.recv()
                    # to retrieve work.
                    # If the Primary does not keep send_request in scope, the send_request buffer can be overwitten by a random process by the time
                    # the replica gets to it.
                    # When the primary does a non-blocking call to retrieve an response form the replica in MPI.Comm.irecv,
                    # it gets another MPI.Request object (rcv_request) whose buffer will contain the response from the replica when it's finished.
                    # The buffers within send_request and rcv_request at separate mem addresses.
                    busy_replica_2_request[replica_rank] = WindowReplicaInfo(replica_rank=replica_rank,
                                                                       work_arguments=window_args,
                                                                       mpi_send_request=send_request,
                                                                       mpi_rcv_request=rcv_request)

                    start_window_nucpos += window_slide

                # Check on replicas
                if busy_replica_2_request:
                    requests = [window_replica_info.mpi_request for window_replica_info in busy_replica_2_request.values()]
                    mpi_status = MPI.Status()
                    idx, err_msg = MPI.Request.waitany(requests=requests, status=mpi_status)
                    done_replica_rank = mpi_status.Get_source()
                    available_replicas.extend([done_replica_rank])
                    busy_replica_2_request.pop(done_replica_rank)
                    if err_msg:
                        LOGGER.error("Received error from replica=" + str(done_replica_rank) + " err_msg=" + str(err_msg))
                    else:
                        LOGGER.debug("Received success from replica=" + str(done_replica_rank))

            LOGGER.debug("Done Launching " + str(total_windows) + " total windows")

            LOGGER.debug("Terminating replicas...")
            for replica_rank in range(1, pool_size):
                comm.isend(obj=None, dest=replica_rank, tag=TAG_TERMINATE)
            LOGGER.debug("Done terminating replicas.")

            LOGGER.debug("About to tabulate results")
            tabulate_results(ref, ref_len, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n,
                             start_nucpos, end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff,
                             output_csv_filename, mode,window_slide, smooth_dist)

            LOGGER.debug("Done tabulating results")

        else:  # replica process does the work
            is_terminated = False
            while not is_terminated:
                try:
                    mpi_status = MPI.Status()
                    # block till the primary tells me to do something
                    window_args = comm.recv(source=PRIMARY_RANK, tag=MPI.ANY_TAG, status=mpi_status)


                    if mpi_status.Get_tag() == TAG_TERMINATE:
                        is_terminated = True
                        LOGGER.debug("Replica of rank %d directed to terminate by primary" % rank)
                    else:
                        str_window_args = ', '.join('{}:{}'.format(key, val) for key, val in window_args.items())
                        LOGGER.debug("Received window_args=" + str_window_args)
                        eval_window(**window_args)
                        comm.send(obj=None, dest=PRIMARY_RANK, tag=TAG_WORK)
                except Exception, e:
                    LOGGER.exception("Failure in replica=" + str(rank))
                    err_msg = traceback.format_exc()
                    comm.send(obj=err_msg, dest=PRIMARY_RANK, tag=TAG_WORK)
    except Exception, e:
        LOGGER.exception("Uncaught Exception.  Aborting")
        comm.Abort()


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("sam_filename", help="full filepath to sam alignment file")
    parser.add_argument("ref", help="Name of reference contig.  Must be one of refnames in SAM file.")
    #parser.add_argument("--ref_len", type=int, help="length of reference contig in nucleotide bases")
    parser.add_argument("--out_dir", default='.',
                        help="output directory in which the pipeline will write all its intermediate files")
    parser.add_argument("--map_qual_cutoff", type=int, default=20,
                        help="mapping quality threshold below which alignments are ignored")
    parser.add_argument("--read_qual_cutoff", type=int, default=20,
                        help="read quality threshold below which bases are converted to Ns")
    parser.add_argument("--max_prop_n", type=float, default=0.1,
                        help="maximum fraction of Ns allowed in the merged paired-end read below which the paired-end"
                             " read is ignored")
    parser.add_argument("--window_size", type=int, default=400, help="window size in nucleotides")
    parser.add_argument("--window_slide", type=int, default=50, help="Number of bases to slide each window by")
    parser.add_argument("--window_breadth_cutoff", type=float, default=0.8,
                        help="fraction of window that merged paired-end read must cover with non-gap and non-N"
                             " nucleotides.  Below this threshold, the read is omitted from the window.")
    parser.add_argument("--window_depth_cutoff", type=int, default=50,
                        help="Minimum number of reads within a valid window.")
    parser.add_argument("--start_nucpos", type=int, default=1,
                        help="1-based start nucleotide position in the reference contig.  The first window will start"
                             " at this position.")
    parser.add_argument("--end_nucpos", type=int, default=None,
                        help="1-based end nucleotide position in the reference contig.  The last window will start at"
                             " or before this position.")
    parser.add_argument("--pvalue", type=float, default=0.05,
                        help=" p-value threshold for determining selection significantly different from neutral"
                             " evolution.")
    parser.add_argument("--threads_per_window", type=int, default=1,
                        help="threads allotted per window.")
    parser.add_argument("--concurrent_windows", type=int, default=1,
                        help="Max number of windows to process concurrently. Ignored when --mpi is defined.")
    parser.add_argument("--output_csv_filename", default='dnds.tsv',
                        help="In DNDS mode, the full filepath of final tab-separated values file containing selection information for"
                             " each codon site in the reference from averaged over multiple windows.")
    parser.add_argument("--hyphy_exe", help="full filepath of HYPHYMP executable.  Default: taken from PATH")
    parser.add_argument("--hyphy_basedir",
                        help="full filepath of HyPhy base directory containing template batch files.  Default:"
                             " /usr/local/lib/hyphy/TemplateBatchFiles/")
    parser.add_argument("--fastree_exe", help="full filepath of FastTreeMP executable.  Default: taken from PATH")
    parser.add_argument("--mode", default='DNDS', help="DNDS: Execute dN/dS analysis for positive (diversifying "
                                                       "selection in codon alignment.  GTR_RATE: Profile "
                                                       "nucleotide substitution rate biases under generalized "
                                                       "non-reversible (6-paramter) model.")

    parser.add_argument("--mpi", action='store_true',
                        help="Runs in MPI mode with multiple processes on multiple nodes. "
                             "If python module mpi4py is not installed, then runs multiple processes on single "
                             "node. Default: false")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    # deep copy of arguments excluding empty values
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
