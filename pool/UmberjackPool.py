"""
Defines worker class that should be subclassed
"""

import logging
import os
import config.settings
import sam.sam_handler as sam_handler
import math
import fasttree.fasttree_handler as fasttree
import hyphy.hyphy_handler as hyphy
import rtt.rtt as rtt
import slice_miseq
import plot.plotter as plotter


config.settings.setup_logging()
LOGGER = logging.getLogger(__name__)

MODE_DNDS = "DNDS"
MODE_GTR_CMP = "GTR_CMP"
MODE_GTR_RATE = "GTR_RATE"
MODE_COUNT_SUBS = "COUNT_SUBS"



def eval_window(sam_filename, ref, out_dir, window_depth_cutoff, window_breadth_cutoff, start_window_nucpos,
            end_window_nucpos, map_qual_cutoff, read_qual_cutoff, max_prop_N, insert, mask_stop_codon, remove_duplicates,
            threads_per_window=config.settings.DEFAULT_THREADS_PER_WINDOW, mode=config.settings.DEFAULT_MODE,
            hyphy_exe=config.settings.DEFAULT_HYPHY_EXE, hyphy_basedir=config.settings.DEFAULT_HYPHY_BASEDIR,
            hyphy_libdir=config.settings.DEFAULT_HYPHY_LIBDIR,
            fastree_exe=config.settings.DEFAULT_FASTTREEMP_EXE,
            debug=False):
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
    :param str hyphy_basedir:  full filepath to directory containing the custom HyPhy template batch files
    :param str hyphy_libdir:  full filepath to directory containing the out of box HyPhy template batch files
    :param str fastree_exe: full filepath to FastTreeMP or FastTree executable
    """
    LOGGER.debug("Eval window {}-{}".format(start_window_nucpos, end_window_nucpos) + " for sam=" + sam_filename + " ref=" + ref)

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
        LOGGER.warn("MSA Window " + msa_window_fasta_filename + " does not satisfy window depth constraints  with " + str(total_slice_seq) + " seq")
    else:
        LOGGER.debug("MSA Window " + msa_window_fasta_filename + " satisfies window depth constraints with " + str(total_slice_seq) + " seq")

        fastree_treefilename = fasttree.make_tree(fasta_fname=msa_window_fasta_filename, threads=threads_per_window, fastree_exe=fastree_exe)

        rooted_treefile = None
        if mode == MODE_COUNT_SUBS:
            rooted_treefile = rtt.make_rooted_tree(unrooted_treefile=fastree_treefilename, threads=threads_per_window)

        if mode == MODE_DNDS:
            hyphy.calc_dnds(codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename,
                            hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, hyphy_libdir=hyphy_libdir,
                            threads=threads_per_window, debug=debug)
        elif mode == MODE_COUNT_SUBS:
            hyphy.count_site_branch_subs(codon_fasta_filename=msa_window_fasta_filename, rooted_treefile=rooted_treefile,
                            hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, hyphy_libdir=hyphy_libdir,
                            threads=threads_per_window, debug=debug)
        elif mode == MODE_GTR_CMP:
            hyphy.calc_nuc_subst(hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, threads=threads_per_window,
                                 codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename)
        elif mode != MODE_GTR_RATE:
            raise  ValueError("Invalid mode " + mode)

    LOGGER.debug("Done Eval window {}-{}".format(start_window_nucpos, end_window_nucpos) + " for sam=" + sam_filename + " ref=" + ref)


class UmberjackPool(object):
    """
    Does the pipeline work.
    Let's subclasses decide how to spread out work amongst child processes.
    """



    def __init__(self, **kwargs):
        self.out_dir = config.settings.DEFAULT_OUT_DIR
        self.map_qual_cutoff = config.settings.DEFAULT_MAP_QUAL_CUTOFF
        self.read_qual_cutoff = config.settings.DEFAULT_READ_QUAL_CUTOFF
        self.max_prop_n = config.settings.DEFAULT_MAX_PROP_N
        self.start_nucpos  = 0
        self.end_nucpos = 0
        self.window_size = 300
        self.window_depth_cutoff = 10
        self.window_breadth_cutoff = 0.875
        self.window_slide = 30
        self.insert = False
        self.mask_stop_codon = True
        self.remove_duplicates = True
        self.output_csv_filename = "./output.csv"
        self.mode = MODE_DNDS
        self.concurrent_windows = 1
        self.threads_per_window = 1
        self.hyphy_exe = config.settings.DEFAULT_HYPHY_EXE
        self.hyphy_basedir = config.settings.DEFAULT_HYPHY_BASEDIR
        self.hyphy_libdir = config.settings.DEFAULT_HYPHY_LIBDIR
        self.fastree_exe = config.settings.DEFAULT_FASTTREEMP_EXE
        self.debug = True
        self.ref = None
        self.sam_filename = None
        self.sam_filename_list = None
        self.mpi = False
        self.__dict__.update(kwargs)


        # if self.mpi:
        #     from pool.MPIPool import  MPIPool
        #     self. = MPIPool(size=None)
        # else:
        #     from pool.SingleNodePool import  SingleNodePool
        #     self. = SingleNodePool(size=self.concurrent_windows)


    def is_parent(self):
        raise NotImplementedError("This object should be subclassed")


    def setup(self):
        """
        Stuff that only the parent should do before child processes work.
        :return
        """
        self.check_sams()
        self.make_sam_ref_outdir()




    def make_sam_ref_outdir(self):
        """
        Makes the output dir or outputdir/sam/ref nested dirs.
        In order to avoid race conditions, make these output directories before spread_working child processes.
        """
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        if self.sam_filename_list or not self.ref:  # multiple sam, ref combos
            for samfile, ref in self.sam_ref_iter():
                sam_ref_outdir = UmberjackPool.get_sam_ref_outdir(pardir=self.out_dir, samfilename=samfile, ref=ref, samfilelist=self.sam_filename_list)

                if not os.path.exists(sam_ref_outdir):
                    os.makedirs(sam_ref_outdir)


    @staticmethod
    def check_sam(sam_filename, ref):
        """
        Checks that either sam_filename is defined and points to valid SAM file(s).
        Checks that all sam file has header and is queryname sorted.
        Checks that if ref isn't specified, then there are refs in the sam header.

        :return bool:  True if everything passes
        :raises ValueError:  if there are issues
        """
        if not sam_handler.is_query_sort(sam_filename):
            LOGGER.fatal("SAM file " + sam_filename + " is not queryname sorted")
            raise ValueError("SAM file " + sam_filename + " is not queryname sorted")

        refs = sam_handler.get_refs(sam_filename)
        if not refs:
            LOGGER.fatal("SAM file " +  sam_filename + " does not have references in header")
            raise ValueError("SAM file " +  sam_filename + " does not have references in header")
        elif ref and ref not in refs:
            LOGGER.fatal("SAM file " +  sam_filename + " does not contain " + ref + " in its references")
            raise ValueError("SAM file " +  sam_filename + " does not contain " + ref + " in its references")

        return True


    def check_sams(self):
        """
        Checks that either sam_filename_list or sam_filename specified
        Checks that all sam files specified have headers and are queryname sorted.
        Checks that if ref isn't specified, then ref is stored in each specified sam file.

        :return bool:  True if everything passes
        :raises ValueError:  if there are issues
        """
        if not self.sam_filename and not self.sam_filename_list:
            LOGGER.fatal("Must specify either sam_filename or sam_filename_list")
            raise ValueError("Must specify either sam_filename or sam_filename_list")
        elif self.sam_filename:
            UmberjackPool.check_sam(self.sam_filename, self.ref)
        else:
            for sam in UmberjackPool.list_from_file(self.sam_filename_list):
                UmberjackPool.check_sam(sam, self.ref)


        return True

    @staticmethod
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
                                                  ref_len=0, do_remove_dup=False)

            LOGGER.debug("Done Full MSA-Fasta from SAM for ref " + ref)
        else:
            LOGGER.warn("Found existing Full MSA-Fasta from SAM for ref " + ref + ".  Not regenerating")

        return msa_fasta_filename

    @staticmethod
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

    @staticmethod
    def list_from_file(listfile):
        """
        Gets a list from newline separated file of items
        :return list[str]:
        """
        items = []
        with open(listfile, 'r') as fh_in:
            for line in fh_in:
                line = line.rstrip()
                items.extend([line])
        return items

    def sam_ref_iter(self):
        """
        Yields an iterator for each sample-ref combination
        :return iter(str, str): Iterator for (samfile, ref) (samfiletuple
        """
        if self.sam_filename:
            queue_samfiles = [self.sam_filename]
        else:
            queue_samfiles = UmberjackPool.list_from_file(self.sam_filename_list)

        for queue_next_samfile in queue_samfiles:
            if self.ref:
                queue_refs = [self.ref]
            else:
                queue_refs = sam_handler.get_refs(queue_next_samfile)

            for queue_next_ref in queue_refs:
                yield queue_next_samfile, queue_next_ref

    @staticmethod
    def get_sam_ref_outdir(pardir, samfilename, ref, samfilelist=None):
        # TODO:  hack if there is only 1 sam and ref, then just spit out the pardir
        if not samfilelist and ref:
            return pardir
        sam_basename = os.path.basename(samfilename)
        sam_prefix = os.path.splitext(sam_basename)[0]
        sam_ref_outdir = pardir + os.sep + sam_prefix + os.sep + ref

        return sam_ref_outdir

    def window_iter(self):
        """
        Yields an iterator for each samfile-ref-window.

        Creates the full MSA fasta and duplicates file for every sam-file-ref.
        :return:
        """
        if self.sam_filename:
            queue_samfiles = [self.sam_filename]
        else:
            queue_samfiles = UmberjackPool.list_from_file(self.sam_filename_list)

        for queue_next_samfile in queue_samfiles:

            if self.ref:
                queue_refs = [self.ref]
            else:
                queue_refs = sam_handler.get_refs(queue_next_samfile)

            for queue_next_ref in queue_refs:
                if not self.start_nucpos:
                    self.start_nucpos = 1
                if not self.end_nucpos:
                    self.end_nucpos = sam_handler.get_reflen(queue_next_samfile, queue_next_ref)

                sam_ref_outdir = UmberjackPool.get_sam_ref_outdir(pardir=self.out_dir, samfilename=queue_next_samfile, ref=queue_next_ref, samfilelist=self.sam_filename_list)

                if self.debug:
                    UmberjackPool.create_full_msa_fasta(sam_filename=queue_next_samfile, out_dir=sam_ref_outdir, ref=queue_next_ref,
                                          mapping_cutoff=self.map_qual_cutoff, read_qual_cutoff=self.read_qual_cutoff,
                                          is_insert=self.insert, is_mask_stop_codon=self.mask_stop_codon)

                # TODO:  only  make this file if debug???
                if self.remove_duplicates:
                    UmberjackPool.create_dup_tsv(sam_filename=queue_next_samfile, out_dir=sam_ref_outdir, ref=queue_next_ref,
                                   mapping_cutoff=self.map_qual_cutoff, read_qual_cutoff=self.read_qual_cutoff,
                                   is_insert=self.insert)

                # All nucleotide positions are 1-based
                total_windows = int(math.ceil(((self.end_nucpos-self.window_size-1) - self.start_nucpos + 1)/self.window_slide))
                LOGGER.debug("Queuing " + str(total_windows) + " total windows for sam=" + queue_next_samfile + " ref=" + queue_next_ref)


                if not self.start_nucpos:
                    queue_ref_start_nucpos = 1
                else:
                    queue_ref_start_nucpos = self.start_nucpos
                if not self.end_nucpos:
                    queue_ref_end_nucpos = sam_handler.get_reflen(queue_next_samfile, queue_next_ref)
                else:
                    queue_ref_end_nucpos =  self.end_nucpos

                for start_window_nucpos in range(queue_ref_start_nucpos, queue_ref_end_nucpos-self.window_size-1, self.window_slide):
                    end_window_nucpos = start_window_nucpos + self.window_size - 1
                    LOGGER.debug("start_window_nucpos=" + str(start_window_nucpos))
                    window_args = {"window_depth_cutoff": self.window_depth_cutoff,
                                   "window_breadth_cutoff": self.window_breadth_cutoff,
                                   "start_window_nucpos": start_window_nucpos,
                                   "end_window_nucpos": end_window_nucpos,
                                   "ref": queue_next_ref,
                                   "out_dir": sam_ref_outdir,
                                   "sam_filename": queue_next_samfile,
                                   "map_qual_cutoff": self.map_qual_cutoff,
                                   "read_qual_cutoff": self.read_qual_cutoff,
                                   "max_prop_N": self.max_prop_n,
                                   "insert": self.insert,
                                   "mask_stop_codon": self.mask_stop_codon,
                                   "remove_duplicates": self.remove_duplicates,
                                   "threads_per_window": self.threads_per_window,
                                   "mode": self.mode,
                                   "hyphy_exe": self.hyphy_exe,
                                   "hyphy_basedir": self.hyphy_basedir,
                                   "hyphy_libdir": self.hyphy_libdir,
                                   "fastree_exe": self.fastree_exe,
                                   "debug": self.debug}
                    yield window_args





    def tabulate_results(self):
        """
        Tabulates results for a samfile & ref into a single csv file.
        Plots results if applicable.
        """
        # TODO:  specify the run parameters in comments
        LOGGER.debug("About to tabulate results")
        if self.mode == MODE_DNDS:
            kwds_list = []
            for samfile, ref in self.sam_ref_iter():
                # TODO:  unhackl
                ref_len = sam_handler.get_reflen(sam_filename=samfile, ref=ref)
                sam_ref_outdir = UmberjackPool.get_sam_ref_outdir(pardir=self.out_dir, samfilename=samfile, ref=ref, samfilelist=self.sam_filename_list)
                if sam_ref_outdir == self.out_dir:
                    output_csv_filename = self.output_csv_filename
                else:
                    output_csv_filename = sam_ref_outdir + os.sep +os.path.basename(samfile).replace(".sam", ".dnds.csv")
                # TODO:  hack to check if m
                kwds_list.append({"dnds_tsv_dir": sam_ref_outdir,
                                  "ref": ref,
                                  "ref_nuc_len": ref_len,
                                  "output_csv_filename": output_csv_filename
                })
            self.spread_work(func=slice_miseq.tabulate_dnds, work_args_iter=kwds_list)

            plot_kws_list = []
            for samfile, ref in self.sam_ref_iter():
                # TODO:  unhackl
                sam_ref_outdir = UmberjackPool.get_sam_ref_outdir(pardir=self.out_dir, samfilename=samfile, ref=ref, samfilelist=self.sam_filename_list)
                if sam_ref_outdir == self.out_dir:
                    output_csv_filename = self.output_csv_filename
                else:
                    output_csv_filename = sam_ref_outdir + os.sep +os.path.basename(samfile).replace(".sam", ".dnds.csv")
                plot_kws_list.append({"dnds_csv": output_csv_filename})

            self.spread_work(func=plotter.plot_dnds, work_args_iter=plot_kws_list)

        elif self.mode == MODE_GTR_RATE:
            kwds_list = []
            for samfile, ref in self.sam_ref_iter():
                # TODO:  unhackl
                sam_ref_outdir = UmberjackPool.get_sam_ref_outdir(pardir=self.out_dir, samfilename=samfile, ref=ref, samfilelist=self.sam_filename_list)
                if sam_ref_outdir == self.out_dir:
                    output_csv_filename = self.output_csv_filename
                else:
                    output_csv_filename = sam_ref_outdir + os.sep +os.path.basename(samfile).replace(".sam", ".gtr.csv")
                kwds_list.append({"fasttree_output_dir": sam_ref_outdir,
                                  "output_csv_filename": output_csv_filename
                })
            self.spread_work(func=slice_miseq.tabulate_rates, work_args_iter=kwds_list)

        LOGGER.debug("Done tabulating results")


    def eval_windows(self):
        """
        spread_work window processes asynchronously.
        :return:
        """
        LOGGER.debug("About to spread window work")
        self.spread_work(func=eval_window, work_args_iter=self.window_iter())
        LOGGER.debug("Done spread window work")

    def start(self):
        """
        Main entry method for the pipeline.
        :return:
        """

        if self.is_parent():
            self.setup()

        # spread_work windows for every sam, ref into the same pool
        self.eval_windows()


        # Tabulate results for every sam, ref
        self.tabulate_results()



    def spread_work(self, func, work_args_iter):
        raise NotImplementedError("This object should be subclassed")