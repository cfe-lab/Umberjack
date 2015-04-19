#!/usr/bin/python
"""
Main entry point for Umberjack.
"""
import os
import sys
import logging
import argparse
import traceback
import math

from sam import sam_handler
import slice_miseq
import pool_traceback
from config_arg_parse import ConfigArgParser
import hyphy.hyphy_handler as hyphy
import fasttree.fasttree_handler as fasttree
import config.settings as settings
import plot.plotter as plotter
from subprocess import CalledProcessError

LOGGER = logging.getLogger(__name__)


# For MPI
PRIMARY_RANK = 0
TAG_WORK = 1
TAG_TERMINATE = 2

MODE_DNDS = "DNDS"
MODE_GTR_CMP = "GTR_CMP"
MODE_GTR_RATE = "GTR_RATE"



def create_full_msa_fasta(sam_filename, out_dir, ref, mapping_cutoff, read_qual_cutoff,
                          is_insert, is_mask_stop_codon):
    """
    Creates a pseudo multiple-sequence aligned fasta file for all reads using pairwise alignment from a SAM file.
    Does not filter based on breadth thresholds or N's.  But it does mask low quality bases and conflicts.

    :return: path to multiple sequence aligned fasta file of all reads
    :rtype : str
    :param str sam_filename: filepath to sam file of read alignments to the reference.  Must have header and must be queryname sorted.
    :param str out_dir: output directory
    :param str ref: reference name
    :param int mapping_cutoff: minimum mapping quality cutoff.  Reads aligned with map quality below this are thrown out.
    :param int read_qual_cutoff: read quality cutoff.  Bases below this cutoff are converted to N's.
    :param bool is_insert:  If True, then keeps insertions with respect to reference
    :param bool is_mask_stop_codon:  If True, then masks stop codons
    """

    sam_filename_nopath = os.path.split(sam_filename)[1]
    sam_filename_prefix = os.path.splitext(sam_filename_nopath)[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    msa_fasta_filename = out_dir + os.sep + sam_filename_prefix + "." + ref + ".msa.fasta"

    LOGGER.debug("Start Full MSA-Fasta from SAM for ref " + ref)
    if not os.path.exists(msa_fasta_filename) or os.path.getsize(msa_fasta_filename) <= 0:
        sam_handler.create_msa_slice_from_sam(sam_filename=sam_filename, ref=ref, out_fasta_filename=msa_fasta_filename,
                                              mapping_cutoff=mapping_cutoff, read_qual_cutoff=read_qual_cutoff,
                                              max_prop_N=1.0, breadth_thresh=0.0, start_pos=0, end_pos=0,
                                              do_insert_wrt_ref=is_insert, do_mask_stop_codon=is_mask_stop_codon,
                                              ref_len=0)

        LOGGER.debug("Done Full MSA-Fasta from SAM for ref " + ref)
    else:
        LOGGER.warn("Found existing Full MSA-Fasta from SAM for ref " + ref + ".  Not regenerating")

    return msa_fasta_filename


def create_dup_tsv(sam_filename, out_dir, ref, mapping_cutoff, read_qual_cutoff, is_insert):
    """
    Spits the read duplicates to tab separated file

    :return: path to multiple sequence aligned fasta file of all reads
    :rtype : str
    :param str sam_filename: filepath to sam file of read alignments to the reference.  Must have header and must be queryname sorted.
    :param str out_dir: output directory
    :param str ref: reference name
    :param int mapping_cutoff: minimum mapping quality cutoff.  Reads aligned with map quality below this are thrown out.
    :param int read_qual_cutoff: read quality cutoff.  Bases below this cutoff are converted to N's.
    :param bool is_insert:  If True, then keeps insertions with respect to reference
    """

    sam_filename_nopath = os.path.split(sam_filename)[1]
    sam_filename_prefix = os.path.splitext(sam_filename_nopath)[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    dup_tsv_filename = out_dir + os.sep + sam_filename_prefix + "." + ref + ".dup.tsv"

    LOGGER.debug("Start Duplicate Reads TSV " + dup_tsv_filename + " from SAM for ref " + ref)
    if not os.path.exists(dup_tsv_filename) or os.path.getsize(dup_tsv_filename) <= 0:
        total_uniq = sam_handler.write_dup_record_tsv(sam_filename=sam_filename, ref= ref,
                                         mapping_cutoff=mapping_cutoff, read_qual_cutoff=read_qual_cutoff, is_insert=is_insert,
                                         out_tsv_filename=dup_tsv_filename)

        LOGGER.debug("Total unique sequences =" + str(total_uniq) + " for " + dup_tsv_filename)
        LOGGER.debug("Done Duplicate Reads TSV for ref " + ref)
    else:
        LOGGER.warn("Found existing Duplicate Reads TSV " + dup_tsv_filename + " for ref " + ref + ".  Not regenerating")

    return dup_tsv_filename


def eval_window(sam_filename, ref, out_dir, window_depth_cutoff, window_breadth_cutoff, start_window_nucpos,
                end_window_nucpos, map_qual_cutoff, read_qual_cutoff, max_prop_N, insert, mask_stop_codon, remove_duplicates,
                threads_per_window=settings.DEFAULT_THREADS_PER_WINDOW, mode=settings.DEFAULT_MODE,
                hyphy_exe=settings.DEFAULT_HYPHY_EXE, hyphy_basedir=settings.DEFAULT_HYPHY_BASEDIR,
                fastree_exe=settings.DEFAULT_FASTTREEMP_EXE):
    """
    Handles the processing for a single window along the genome.
    Creates the multiple sequence aligned fasta file for the window.
    Feeds the window multiple-sequence aligned fasta file to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.

    :param str sam_filename: full file path to sam file
    :param str ref: reference name
    :param str out_dir: full path to directory that will hold the results and intermediate files
    :param int window_depth_cutoff:  the minimum number of required reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_cutoff: the minimum fraction of a window that merged paired-end read must cover to be included in the window.
    :param int start_window_nucpos:  1-based start nucleotide position of the window
    :param int end_window_nucpos:  1-based end nucleotide position of the window
    :param int map_qual_cutoff:  mapping quality threshold below which read is thrown out
    :param int read_qual_cutoff:  read base quality threshold below which base is masked with N
    :param float max_prop_N:  maximum fraction of single-end read or merged paired read above which read is thrown out
    :param bool insert:  if True, inserts with respect to reference kept
    :param bool mask_stop_codon:  if True, stop codons masked
    :param bool remove_duplicates:  If True, then removes duplicate sequence from the window multiple sequence alignment fasta.
        Sequences are only considered duplicates if they start on the same coordinate with respect to the reference as another sequence,
        and both sequences have matching bases, gaps, N's after quality  masking, insertion processing but before stop codon masking.
    :param int threads_per_window: number of threads allotted to processing this window  (only FastTree and HyPhy will be multithreaded)
    :param str mode: one of [DNDS, GTR_RATE].  Calculate codon site dN/dS or window-specific nucleotide substitution rates.
    :param str hyphy_exe: full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe: full filepath to FastTreeMP or FastTree executable
    """

    LOGGER.debug("sam_filename=" + str(sam_filename) + "\n" +
                 "ref=" + ref +
                 "start_nucpos=" + str(start_window_nucpos) + "\n" +
                 "end_nucpos=" + str(end_window_nucpos) )

    sam_filename_nopath = os.path.split(sam_filename)[1]
    sam_filename_prefix = os.path.splitext(sam_filename_nopath)[0]

    msa_window_filename_prefix = out_dir + os.sep + sam_filename_prefix + "." + str(start_window_nucpos) + "_" + str(end_window_nucpos)
    msa_window_fasta_filename = msa_window_filename_prefix + ".fasta"
    total_slice_seq = sam_handler.create_msa_slice_from_sam(sam_filename=sam_filename, ref=ref,
                                                            out_fasta_filename=msa_window_fasta_filename,
                                                            mapping_cutoff=map_qual_cutoff,
                                                            read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N,
                                                            breadth_thresh=window_breadth_cutoff,
                                                            start_pos=start_window_nucpos, end_pos=end_window_nucpos,
                                                            do_insert_wrt_ref=insert,
                                                            do_mask_stop_codon=mask_stop_codon,
                                                            do_remove_dup=remove_duplicates)


    # Check whether the msa sliced fasta has enough reads to make a good tree
    if total_slice_seq < window_depth_cutoff:
        LOGGER.warn("MSA Window " + msa_window_fasta_filename + " does not satisfy window depth constraints")
    else:
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


def tabulate_results(ref, sam_filename, out_dir, output_csv_filename, mode, **kwargs):
    """
    Tabulates results into a single csv file
    :param str ref:  reference name
    :param str sam_filename:  full file path to sam file
    :param str out_dir:  full file path to output directory in which results csv will be written.
    :param str output_csv_filename:  name of results csv file
    :param str mode:  one of [DNDS, GTR_RATE] for calculating per-codon site dN/dS  averaged over overlapping windows or
                    calculating window specific general time reversible nucleotide substitution rates.
    :param kwargs:  parameters used to generate the results.  Will be written as comments in the results csv file.
    """
    ref_len = sam_handler.get_reflen(sam_filename, ref)
    comments = ("ref=" + ref + "," +
                "ref_len=" + str(ref_len) + "," +
                "sam=" + sam_filename + ","
                "out_dir=" + out_dir + ","
                "mode=" + mode + ",")
    for key, value in kwargs.iteritems():
        comments +=  "{}={}, ".format(key, value)

    if mode == MODE_DNDS:
        LOGGER.debug("Start Ave Dn/DS for all windows for ref " + ref + " " + output_csv_filename)
        seq_dnds_info = slice_miseq.tabulate_dnds(dnds_tsv_dir=out_dir, ref=ref, ref_nuc_len=ref_len,
                                                  output_csv_filename=output_csv_filename, comments=comments)
        LOGGER.debug("Done Ave Dn/DS for all windows  for ref " + ref + ".  Wrote to " + output_csv_filename)
        return seq_dnds_info
    elif mode == MODE_GTR_RATE:
        LOGGER.debug("Start Tabulate GTR Rates for All Windows For Ref " + ref + " " + output_csv_filename)
        slice_miseq.tabulate_rates(fasttree_output_dir=out_dir, output_csv_filename=output_csv_filename, comments=comments)
        LOGGER.debug("Done Tabulate GTR Rates for all windows for ref " + ref + " " + output_csv_filename)
    else:
        LOGGER.debug("Done all windows  for ref " + ref)


def plot_results(output_csv, mode=MODE_DNDS):
    """
    Plot average of site-wise results along genome.
    :param str output_csv:  path to output csv file
    :param str mode:  Only handles DNDS right now
    :return:
    """
    try:
        if mode == MODE_DNDS:
            LOGGER.debug("Plotting Ave dN/dS results for " + output_csv)
            plotter.plot_dnds(output_csv)
        else:
            LOGGER.warn("No plots available for mode " + mode)
    except CalledProcessError, e:
        LOGGER.exception("Unable to plot dnds values to pdf. \n", e.message)


def eval_windows_async(ref, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                       end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, window_slide,
                       insert, mask_stop_codon, remove_duplicates,
                       output_csv_filename, mode=settings.DEFAULT_MODE,
                       threads_per_window=1, concurrent_windows=1,
                       hyphy_exe=settings.DEFAULT_HYPHY_EXE, hyphy_basedir=settings.DEFAULT_HYPHY_BASEDIR, fastree_exe=settings.DEFAULT_FASTTREEMP_EXE,
                       debug=False):
    """
    Launch a separate process to analyze each window.
    Each window can use up to <threads_per_window> threads.

    The maximum number of CPU cores required on a single node = threads_per_window x concurrent_windows.
    :param str ref:  reference name
    :param str sam_filename:  full filepath to sam file of read alignments against the reference
    :param str out_dir:  output directory
    :param int map_qual_cutoff:  mapping quality threshold below which alignments are thrown out
    :param int read_qual_cutoff:  read quality threshold below which bases are converted to N's
    :param float max_prop_n:  maximum fraction of merged paired-end read that are N's within the slice.  Below this threshold the read pair is thrown out.
    :param int start_nucpos:  1-based start nucleotide position with respect to the reference for the window
    :param int end_nucpos:  1-based end nucleotide position with respect to the reference for the window
    :param int window_size:  size of each window in nucleotide bases
    :param int window_depth_cutoff:  the minimum number of required single-end reads or merged paired reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_cutoff:  the minimum fraction of a window a single-end read or merged paired read must cover to be included in the window.
    :param int window_slide:  Number of nucleotide bases to slide the window by.  Default=30.
    :param bool insert:  if True, then keeps insertions with respect to the reference.
    :param bool mask_stop_codon:  If True, then masks stop codons in the window multiple sequence alignment fasta
    :param bool remove_duplicates:  If True, then removes duplicate sequence from the window multiple sequence alignment fasta.
        Sequences are only considered duplicates if they start on the same coordinate with respect to the reference as another sequence,
        and both sequences have matching bases, gaps, N's after quality  masking, insertion processing but before stop codon masking.
    :param int threads_per_window:  number of threads allotted to processing a single window  (only FastTree and HyPhy will be multithreaded).  Default 1.
    :param int concurrent_windows:  the number of windows to process at the same time.
    :param str output_csv_filename:  name of output comma separated file containing per-site results averaged across overlapping windows .  Will be created under out_dir.
    :param str mode:  one of [DNDS, GTR_RATE] for calculating codon-site specific dN/dS, or window-specific general time reversible nucleotide substitution rates
    :param str hyphy_exe:  full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe:  full filepath to FastTreeMP executable
    :param bool debug:  if True, outputs full genome multiple sequence alignment
    """
    fcn_args = locals()

    if not start_nucpos:
        start_nucpos = 1
    if not end_nucpos:
        end_nucpos = sam_handler.get_reflen(sam_filename, ref)

    pool = pool_traceback.LoggingPool(processes=concurrent_windows)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if debug:
        # Create a pseudo multiple-sequence aligned fasta file
        create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref,
                              mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
                              is_insert=insert, is_mask_stop_codon=mask_stop_codon)

    if remove_duplicates:
        create_dup_tsv(sam_filename=sam_filename, out_dir=out_dir, ref=ref,
                       mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
                       is_insert=insert)

    # All nucleotide positions are 1-based
    total_windows = int(math.ceil(((end_nucpos - window_size + 1) - start_nucpos + 1)/window_slide))
    LOGGER.debug("There are " + str(total_windows) + " total windows to process")


    process_results = []
    for start_window_nucpos in range(start_nucpos, end_nucpos-window_size+2, window_slide):
        end_window_nucpos = min(start_window_nucpos + window_size - 1, end_nucpos)
        window_args = {"window_depth_cutoff": window_depth_cutoff,
                       "window_breadth_cutoff": window_breadth_cutoff,
                       "start_window_nucpos": start_window_nucpos,
                       "end_window_nucpos": end_window_nucpos,
                       "ref": ref,
                       "out_dir": out_dir,
                       "sam_filename": sam_filename,
                       "map_qual_cutoff": map_qual_cutoff,
                       "read_qual_cutoff": read_qual_cutoff,
                       "max_prop_N": max_prop_n,
                       "insert": insert,
                       "mask_stop_codon": mask_stop_codon,
                       "remove_duplicates": remove_duplicates,
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


    LOGGER.debug("About to tabulate results")
    tabulate_results(**fcn_args)
    LOGGER.debug("Done tabulating results")

    LOGGER.debug("About to plot results")
    plot_results(output_csv=output_csv_filename)
    LOGGER.debug("Done plot results")



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


def eval_windows_mpi(ref, sam_filename, out_dir, map_qual_cutoff, read_qual_cutoff, max_prop_n, start_nucpos,
                     end_nucpos, window_size, window_depth_cutoff, window_breadth_cutoff, window_slide,
                     insert, mask_stop_codon, remove_duplicates,
                     output_csv_filename, mode, threads_per_window,
                     hyphy_exe, hyphy_basedir,
                     fastree_exe, debug):
    """
    Launch a separate process to analyze each window via MPI.  Similar to eval_windows_async, but uses MPI.

    :param str ref:  reference name
    :param str sam_filename:  full filepath to sam file of read alignments against the reference
    :param str out_dir:  output directory
    :param int map_qual_cutoff:  mapping quality threshold below which alignments are thrown out
    :param int read_qual_cutoff:  read quality threshold below which bases are converted to N's
    :param float max_prop_n:  maximum fraction of merged paired-end read that are N's within the slice.  Below this threshold the read pair is thrown out.
    :param int start_nucpos:  1-based start nucleotide position with respect to the reference for the window
    :param int end_nucpos:  1-based end nucleotide position with respect to the reference for the window
    :param int window_size:  size of each window in nucleotide bases
    :param int window_depth_cutoff:  the minimum number of required single-end reads or merged paired reads that meet the breadth threshold below which the window is thrown out
    :param float window_breadth_cutoff:  the minimum fraction of a window a single-end read or merged paired read must cover to be included in the window.
    :param int window_slide:  Number of nucleotide bases to slide the window by.
    :param bool insert:  if True, then keeps insertions with respect to the reference.
    :param bool mask_stop_codon:  If True, then masks stop codons in the window multiple sequence alignment fasta.
    :param bool remove_duplicates:  If True, then removes duplicate sequence from the window multiple sequence alignment fasta.
        Sequences are only considered duplicates if they start on the same coordinate with respect to the reference as another sequence,
        and both sequences have matching bases, gaps, N's after quality  masking, insertion processing but before stop codon masking.
    :param int threads_per_window:  number of threads allotted to processing a single window  (only FastTree and HyPhy will be multithreaded).  Default 1.
    :param str output_csv_filename:  name of output comma separated file containing per-site results averaged across overlapping windows .  Will be created under out_dir.
    :param str mode:  one of [DNDS, GTR_RATE] for calculating codon-site specific dN/dS, or window-specific general time reversible nucleotide substitution rates.
    :param str hyphy_exe:  full filepath to HYPHYMP executable
    :param str hyphy_basedir:  full filepath to HyPhy base directory containing the template batch files
    :param str fastree_exe:  full filepath to FastTreeMP executable
    :param bool debug:  if True, outputs full genome multiple sequence alignment
    """
    fcn_args = locals()

    if not start_nucpos:
        start_nucpos = 1
    if not end_nucpos:
        end_nucpos = sam_handler.get_reflen(sam_filename, ref)

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    try:

        rank = comm.Get_rank()
        LOGGER.debug("I am rank=" + str(rank))
        LOGGER.debug("I am on machine=" + str(MPI.Get_processor_name()))

        if rank == PRIMARY_RANK:
            pool_size = comm.Get_size()
            LOGGER.debug("Pool size = " + str(pool_size))

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            if debug:
                create_full_msa_fasta(sam_filename=sam_filename, out_dir=out_dir, ref=ref,
                                      mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
                                      is_insert=insert, is_mask_stop_codon=mask_stop_codon)

            if remove_duplicates:
                create_dup_tsv(sam_filename=sam_filename, out_dir=out_dir, ref=ref,
                       mapping_cutoff=map_qual_cutoff, read_qual_cutoff=read_qual_cutoff,
                       is_insert=insert)

            # All nucleotide positions are 1-based
            total_windows = int(math.ceil(((end_nucpos-window_size+1) - start_nucpos + 1)/window_slide))
            LOGGER.debug("Launching " + str(total_windows) + " total windows")

            available_replicas = range(1, pool_size)
            busy_replica_2_request = {}

            start_window_nucpos = start_nucpos


            # TODO:  instead of polling, use a callback function whenever a replicate returns
            while start_window_nucpos <= end_nucpos - window_size +1 or busy_replica_2_request:
                LOGGER.debug("busy_replica_2_request=" + str(busy_replica_2_request))
                LOGGER.debug("start_window_nucpos=" + str(start_window_nucpos))
                # Assign work to replicas
                while start_window_nucpos <= end_nucpos-window_size+1 and available_replicas:
                    end_window_nucpos = min(start_window_nucpos + window_size - 1, end_nucpos)

                    window_args = {"window_depth_cutoff": window_depth_cutoff,
                                   "window_breadth_cutoff": window_breadth_cutoff,
                                   "start_window_nucpos": start_window_nucpos,
                                   "end_window_nucpos": end_window_nucpos,
                                   "ref": ref,
                                   "out_dir": out_dir,
                                   "sam_filename": sam_filename,
                                   "map_qual_cutoff": map_qual_cutoff,
                                   "read_qual_cutoff": read_qual_cutoff,
                                   "max_prop_N": max_prop_n,
                                   "insert": insert,
                                   "mask_stop_codon": mask_stop_codon,
                                   "remove_duplicates": remove_duplicates,
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
                    LOGGER.debug(str(start_window_nucpos))

                # Check on replicas
                if busy_replica_2_request:
                    requests = [window_replica_info.mpi_rcv_request for window_replica_info in busy_replica_2_request.values()]
                    mpi_status = MPI.Status()
                    idx, err_msg = MPI.Request.waitany(requests, mpi_status)
                    done_replica_rank = mpi_status.Get_source()
                    available_replicas.extend([done_replica_rank])
                    del busy_replica_2_request[done_replica_rank]
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
            tabulate_results(**fcn_args)
            LOGGER.debug("Done tabulating results")

            LOGGER.debug("About to plot results")
            plot_results(output_csv=output_csv_filename)
            LOGGER.debug("Done plot results")

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
                except Exception:
                    LOGGER.exception("Failure in replica=" + str(rank))
                    err_msg = traceback.format_exc()
                    comm.send(obj=err_msg, dest=PRIMARY_RANK, tag=TAG_WORK)
    except Exception:
        LOGGER.exception("Uncaught Exception.  Aborting")
        comm.Abort()


def main():
    """
    Parses commandline arguments and kicks off mpi or multiprocessing versions of window evaluations.
    """
    parser = ConfigArgParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                             description="Version " + str(settings.VERSION),
                             fromfile_prefix_chars='@')
    parser.add_argument("-f", help="Full filepath to config file.  If defined, then ignores all other commandline arguments.")
    parser.add_argument("--sam_filename", help="Full filepath to SAM alignment file.  Must be queryname sorted and must have header.")
    parser.add_argument("--ref", help="Name of reference contig.  Must be one of references in SAM file.")
    parser.add_argument("--out_dir", default='./',
                        help="Output directory in which the pipeline will write all its results and intermediate files")
    parser.add_argument("--map_qual_cutoff", type=int, default=settings.DEFAULT_MAP_QUAL_CUTOFF,
                        help="Mapping quality threshold below which alignments are ignored.")
    parser.add_argument("--read_qual_cutoff", type=int, default=settings.DEFAULT_READ_QUAL_CUTOFF,
                        help="Read quality threshold below which bases are converted to Ns.")
    parser.add_argument("--max_prop_n", type=float, default=settings.DEFAULT_MAX_PROP_N,
                        help="Maximum fraction of Ns allowed in single-end reads or merged paired reads below which the read is ignored.")
    parser.add_argument("--window_size", type=int, default=settings.DEFAULT_WINDOW_SIZE,
                        help="Window size in nucleotides.")
    parser.add_argument("--window_slide", type=int, default=settings.DEFAULT_WINDOW_SLIDE,
                        help="Number of bases to slide each window by.")
    parser.add_argument("--window_breadth_cutoff", type=float, default=settings.DEFAULT_WINDOW_BREADTH_CUTOFF,
                        help="fraction of window that a single-end or merged paired read must cover with non-gap and non-N"
                             " nucleotides.  Below this threshold, the read is omitted from the window.")
    parser.add_argument("--window_depth_cutoff", type=int, default=settings.DEFAULT_WINDOW_DEPTH_CUTOFF,
                        help="Minimum number of reads within a valid window.")
    parser.add_argument("--start_nucpos", type=int, default=1,
                        help="1-based start nucleotide position in the reference contig.  The first window will start"
                             " at this position.  Default: 1")
    parser.add_argument("--end_nucpos", type=int, default=0,
                        help="1-based end nucleotide position in the reference contig.  The last window will start at"
                             " or before this position.  If 0, then automatically set to last position in the reference contig.")
    parser.add_argument("--insert",  action='store_true',
                        help="Whether to keep insertions with respect to the reference.")
    parser.add_argument("--mask_stop_codon",  action='store_true',
                        help="Whether to mask stop codons with NNN.  Automatically set to True when mode is DNDS")
    parser.add_argument("--remove_duplicates",  action='store_true',
                        help="Whether to remove duplicate sequences.  Automatically set to True when mode is DNDS.  "
                             "Sequences are only considered duplicates if they start on the same coordinate " +
                             "with respect to the reference, and both sequences have " +
                             "matching bases, gaps, N's after quality  masking and insertion processing + "
                             "but before stop codon masking.  The entire sequence, not just the portion that fits " +
                             " into the window is compared for duplication.  " +
                             " The sequence with the highest sum of bases (ACGT) aligned to the"
                             " reference will be written to the window fastas.  The names of all duplicated reads will be " +
                             " written to a tab-separated file under <out_dir>")
    parser.add_argument("--threads_per_window", type=int, default=settings.DEFAULT_THREADS_PER_WINDOW,
                        help="threads allotted per window.")
    parser.add_argument("--concurrent_windows", type=int, default=settings.DEFAULT_CONCURRENT_WINDOWS,
                        help="Max number of windows to process concurrently. Ignored when --mpi is defined.")
    parser.add_argument("--output_csv_filename", default=settings.DEFAULT_OUTPUT_CSV_FILENAME,
                        help="In DNDS mode, the full filepath of final comma-separated values file containing selection information for"
                             " each codon site in the reference from averaged over overlapping windows.  " +
                             "In GTR_RATE mode, the full filepath of final comma-separated values file containing " +
                             "window-specific general time reversible nucleotide substitution rates.")
    parser.add_argument("--hyphy_exe", default=settings.DEFAULT_HYPHY_EXE,
                        help="full filepath of HYPHYMP executable.  Default: taken from PATH")
    parser.add_argument("--hyphy_basedir", default=settings.DEFAULT_HYPHY_BASEDIR,
                        help="full filepath of HyPhy base directory containing template batch files.")
    parser.add_argument("--fastree_exe", default=settings.DEFAULT_FASTTREEMP_EXE,
                        help="full filepath of FastTreeMP or FastTree executable.  Default: taken from PATH")
    parser.add_argument("--mode", default=settings.DEFAULT_MODE, choices=[MODE_DNDS, MODE_GTR_RATE],
                        help="DNDS: Execute dN/dS analysis for positive (diversifying "
                             "selection in codon alignment.  GTR_RATE: Profile "
                             "nucleotide substitution rate biases under generalized "
                             "non-reversible (6-parameter) model.")

    parser.add_argument("--mpi", action='store_true',
                        help="Runs in MPI mode with multiple processes on multiple nodes. "
                             "If python module mpi4py is not installed, then runs multiple processes on single "
                             "node.")

    parser.add_argument("--debug", action='store_true',
                        help="Whether to keep all intermediate files, generate full genome multiple sequence alignment, set debug logging."
                             "Overrides the logging level configured in logging.conf.")


    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    if args.f:
        if len(sys.argv) > 2:
            LOGGER.warn("Using config file " + args.f + ".  Ignoring all other commandline arguments.")
        args = parser.parse_args([parser.fromfile_prefix_chars +  args.f])


    # Automatically set mask_stop_codon to True for dN/dS analysis.
    if args.mode == MODE_DNDS:
        args.mask_stop_codon = True
        args.remove_duplicates = True
        LOGGER.warn("Auto-setting --mask_stop_codon to True since Umberjack is in " + MODE_DNDS + " mode")
        LOGGER.warn("Auto-setting --remove_duplicates to True since Umberjack is in " + MODE_DNDS + " mode")

    # deep copy of arguments excluding empty values
    args_msg = ""
    eval_windows_args = {}
    for key, val in vars(args).iteritems():
        args_msg += "{}={} ".format(key, val)
        if val is not None:
            eval_windows_args[key] = val
    LOGGER.info("Launched umberjack.py with settings: " + args_msg)


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



    # Clean out the commandline arguments not used in the eval_windows* methods.
    eval_windows_args.pop("mpi", None)
    eval_windows_args.pop("f", None)

    if do_mpi:
        eval_windows_mpi(**eval_windows_args)
    else:
        eval_windows_async(**eval_windows_args)


if __name__ == '__main__':
    settings.setup_logging()
    main()
