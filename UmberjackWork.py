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
import Utility
from collections import namedtuple
import csv
import sam.sam_handler
import sys

config.settings.setup_logging()
LOGGER = logging.getLogger(__name__)

MODE_DNDS = "DNDS"
MODE_GTR_CMP = "GTR_CMP"
MODE_GTR_RATE = "GTR_RATE"
MODE_COUNT_SUBS = "COUNT_SUBS"




# Function needs to be at level of module to be called via multiprocessing.  Can't be defined in class inside module.
def eval_window(output_prefix, window_depth_cutoff, window_breadth_cutoff, start_window_nucpos, end_window_nucpos,
                map_qual_cutoff, read_qual_cutoff, max_prop_N, insert, mask_stop_codon, remove_duplicates,
                sam_filename=None, ref=None, msa_fasta=None,
                threads_per_window=config.settings.DEFAULT_THREADS_PER_WINDOW, mode=config.settings.DEFAULT_MODE,
                hyphy_exe=config.settings.DEFAULT_HYPHY_EXE, hyphy_basedir=config.settings.DEFAULT_HYPHY_BASEDIR,
                hyphy_libdir=config.settings.DEFAULT_HYPHY_LIBDIR, fastree_exe=config.settings.DEFAULT_FASTTREEMP_EXE,
                debug=False):
    """
    Handles the processing for a single window along the genome.
    Creates the multiple sequence aligned fasta file for the window.
    Feeds the window multiple-sequence aligned fasta file to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.

    :param str sam_filename: full file path to sam file
    :param str ref: reference name
    :param str msa_fasta:  full filepath to multiple sequence aligned fasta file.
            Must specify either [sam_filename & ref] or msa_fasta
    :param str output_prefix: file path and filename prefix of output.  Do not include file suffix or any window specifiers.
        Window specifiers and file suffixes will be appended accordingly.  Any directories for the output should already exist.
        E.G.  ouput_prefix = /MyOutputDir/MyOutputFilePrefix
        will yield intermediate files:
            /MyOutputDir/MyOutputFilePrefix.1_300.fasta
            /MyOutputDir/MyOutputFilePrefix.1_300.nwk
        and will yield final result files:
            /MyOutputDir/MyOutputFilePrefix.dnds.csv

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
    :param bool debug:  whether to use DEBUG logging and output potentially large
        files useful for debugging such as multiple-sequence aligned fasta for entire reference length.
    """
    LOGGER.debug("Eval window {}-{}".format(start_window_nucpos, end_window_nucpos) +
                 " for sam=" + str(sam_filename) + " ref=" + str(ref) + " msa_fasta=" + str(msa_fasta))
    if ((sam_filename and ref and msa_fasta) or
            (not msa_fasta and (not sam_filename or not ref))):
        raise ValueError("Must specify just msa_fasta or (sam_filename and ref)")

    elif sam_filename and ref:
        msa_window_filename_prefix = output_prefix + "." + str(start_window_nucpos) + "_" + str(end_window_nucpos)
        msa_window_fasta_filename = msa_window_filename_prefix + ".fasta"

        total_slice_seq = sam.sam_handler.create_msa_slice_from_sam(sam_filename=sam_filename, ref=ref,
                                                                out_fasta_filename=msa_window_fasta_filename,
                                                                mapping_cutoff=map_qual_cutoff,
                                                                read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N,
                                                                breadth_thresh=window_breadth_cutoff,
                                                                start_pos=start_window_nucpos, end_pos=end_window_nucpos,
                                                                do_insert_wrt_ref=insert,
                                                                do_mask_stop_codon=mask_stop_codon,
                                                                do_remove_dup=remove_duplicates)

    else:
        msa_window_filename_prefix = output_prefix + "." + str(start_window_nucpos) + "_" + str(end_window_nucpos)
        msa_window_fasta_filename = msa_window_filename_prefix + ".fasta"

        # TODO:  handle duplicates
        total_slice_seq = slice_miseq.create_slice_msa_fasta(fasta_filename=msa_fasta,
                                                         out_fasta_filename=msa_window_fasta_filename,
                                                         start_pos=start_window_nucpos,
                                                         end_pos=end_window_nucpos,
                                                         max_prop_N=max_prop_N,
                                                         breadth_thresh=window_breadth_cutoff,
                                                         do_mask_stop_codon=mask_stop_codon)


    # Check whether the msa sliced fasta has enough reads to make a good tree
    if total_slice_seq < window_depth_cutoff:
        LOGGER.warn("MSA Window " + msa_window_fasta_filename + " does not satisfy window depth constraints  with " + str(total_slice_seq) + " seq")
    else:
        LOGGER.debug("MSA Window " + msa_window_fasta_filename + " satisfies window depth constraints with " + str(total_slice_seq) + " seq")

        if mode == MODE_DNDS:
            fastree_treefilename = fasttree.make_tree(fasta_fname=msa_window_fasta_filename, threads=threads_per_window, fastree_exe=fastree_exe)
            hyphy.calc_dnds(codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename,
                            hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, hyphy_libdir=hyphy_libdir,
                            threads=threads_per_window, debug=debug)

        elif mode == MODE_COUNT_SUBS:
            if rtt.is_timeable(msa_window_fasta_filename):
                fastree_treefilename = fasttree.make_tree(fasta_fname=msa_window_fasta_filename, threads=threads_per_window, fastree_exe=fastree_exe)
                rooted_treefile = rtt.make_rooted_tree(unrooted_treefile=fastree_treefilename, threads=threads_per_window)
                hyphy.count_site_branch_subs(codon_fasta_filename=msa_window_fasta_filename, rooted_treefile=rooted_treefile,
                            hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, hyphy_libdir=hyphy_libdir,
                            threads=threads_per_window, debug=debug)
            else:
                LOGGER.warn("MSA window fasta " + msa_window_fasta_filename +
                            " does not contain multiple timepoints.  Unable to root tree or time ancestor.")

        elif mode == MODE_GTR_CMP:
            fastree_treefilename = fasttree.make_tree(fasta_fname=msa_window_fasta_filename, threads=threads_per_window, fastree_exe=fastree_exe)
            hyphy.calc_nuc_subst(hyphy_exe=hyphy_exe, hyphy_basedir=hyphy_basedir, threads=threads_per_window,
                                 codon_fasta_filename=msa_window_fasta_filename, tree_filename=fastree_treefilename)
        elif mode != MODE_GTR_RATE:
            raise  ValueError("Invalid mode " + mode)

    LOGGER.debug("Done Eval window {}-{}".format(start_window_nucpos, end_window_nucpos) +
                 " for sam=" + str(sam_filename) + " ref=" + str(ref) + " msa_fasta=" + str(msa_fasta))



class UmberjackWork(object):
    """
    Does the pipeline work.
    Let's internal process pool class decide how to spread out work amongst child processes.
    """

    Input_CSV_Fieldnames = ["File", "Ref", "OutputPrefix"]
    Job = namedtuple("Job", field_names=Input_CSV_Fieldnames + ["FileType"])
    class InputFileTypes:
        SAM = "SAM"
        MSA = "MSA"

    def __init__(self, **kwargs):
        self.out_dir = config.settings.DEFAULT_OUT_DIR
        self.map_qual_cutoff = config.settings.DEFAULT_MAP_QUAL_CUTOFF
        self.read_qual_cutoff = config.settings.DEFAULT_READ_QUAL_CUTOFF
        self.max_prop_n = config.settings.DEFAULT_MAX_PROP_N
        self.start_nucpos  = 0
        self.end_nucpos = 0
        self.window_size = config.settings.DEFAULT_WINDOW_SIZE
        self.window_depth_cutoff = config.settings.DEFAULT_WINDOW_DEPTH_CUTOFF
        self.window_breadth_cutoff = config.settings.DEFAULT_WINDOW_BREADTH_CUTOFF
        self.window_slide = config.settings.DEFAULT_WINDOW_SLIDE
        self.insert = config.settings.DEFAULT_INSERT
        self.mask_stop_codon = config.settings.DEFAULT_MASK_STOP_CODON
        self.remove_duplicates = config.settings.DEFAULT_REMOVE_DUPLICATES
        self.mode = MODE_DNDS
        self.concurrent_windows = config.settings.DEFAULT_CONCURRENT_WINDOWS
        self.threads_per_window = config.settings.DEFAULT_THREADS_PER_WINDOW
        self.hyphy_exe = config.settings.DEFAULT_HYPHY_EXE
        self.hyphy_basedir = config.settings.DEFAULT_HYPHY_BASEDIR
        self.hyphy_libdir = config.settings.DEFAULT_HYPHY_LIBDIR
        self.fastree_exe = config.settings.DEFAULT_FASTTREEMP_EXE
        self.debug = config.settings.DEFAULT_DEBUG
        self.ref = None
        self.sam_filename = None
        self.msa_fasta = None
        self.input_csv = None
        self.mpi = False
        self.pool = None
        self.__dict__.update(kwargs)


    def setup(self):
        """
        Stuff that only the parent should do before child processes work.
        :return
        """
        if self.pool.is_parent():
            self.check_input()
            self.make_outdir()




    def make_outdir(self):
        """
        Makes all the required output directories.
        In order to avoid race conditions, make these output directories before spreading work to child processes.
        """
        # multiple sam, ref combos
        for samfile, ref in self.sam_ref_iter():
            sam_ref_outdir = UmberjackWork.get_auto_outdir(pardir=self.out_dir, input_file=samfile, ref=ref)
            if not os.path.exists(sam_ref_outdir):
                os.makedirs(sam_ref_outdir)

        if self.msa_fasta:
            outdir = UmberjackWork.get_auto_outprefix(pardir=self.out_dir, input_file=self.msa_fasta, ref=self.ref)
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        # multiple sam or msa fasta, ref combos
        for job in self.input_csv_iter():
            outdir = os.path.dirname(job.OutputPrefix)
            if not os.path.exists(outdir):
                os.makedirs(outdir)



    def check_input_csv(self):
        """
        If self.input_csv is defined, then checks if the Input CSV has
        headers 'File,Ref,OutputPrefix'.
        Checks that files specified are valid sam files or valid multiple-sequence aligned (MSA) fastas.
        :return bool:  true if everything is OK
        :raises ValueError:  if the input csv is missing or empty,
         if any of the entries in the input csv is incomplete,
         if any of the sam or msa fasta files are invalid.
        """
        if not self.input_csv:
            raise ValueError("Undefined input csv")

        if not os.path.exists(self.input_csv) or not os.path.getsize(self.input_csv):
            raise ValueError("Input CSV " + self.input_csv + " does not exist or is empty")

        with open(self.input_csv, 'rU') as fh_in:
            reader = csv.DictReader(fh_in)
            missing_fields = [x for x in UmberjackWork.Input_CSV_Fieldnames if x not in reader.fieldnames]
            if len(missing_fields) > 0:
                raise ValueError("Input CSV " + self.input_csv + " missing fields " + ",".join(missing_fields))

            # 1-based line numbers.  Header of input csv  counts as line 1.
            for i, row in enumerate(reader):
                filename = row["File"]
                ref = row["Ref"]
                output_prefix = row["OutputPrefix"]
                if not filename or not ref or not output_prefix:
                    raise ValueError("Line " + str(i + 2) + " is missing information in " + self.input_csv)

                try:
                    if UmberjackWork.is_sam(filename):
                        UmberjackWork.check_sam(sam_filename=filename, ref=ref)

                    elif UmberjackWork.is_msa_fasta(filename):
                        UmberjackWork.check_msa_fasta(msa_fasta=filename)
                    else:
                        raise ValueError("File " + filename + " is neither sam or multiple sequence aligned fasta format")

                except Exception, e:
                    raise (ValueError("File " + filename + " specified on line " + str(i+2) +
                                      " in " + self.input_csv + " is not valid.\n"),  # exception type
                           None,  # exception value
                           sys.exc_info()[:2]  # previous exception stacktrace
                    )



    def input_csv_window_iter(self):
        """
        Parses self.input_csv to get windows to process.
        Does not do any error checking.  Use UmberjackWork.check_input_csv() for error checking.
        Must have headers 'File,Ref,OutputPrefix'.  Should be comma-delimited.
        :yields iterator of dict: For each (samfile, ref, output prefix) or (msa fasta, ref, output prefix) combination,
            iterates over all the windows and creates dict arguments that can be passed to eval_window()
        """
        if self.input_csv:
            with open(self.input_csv, 'rU') as fh_in:
                reader = csv.DictReader(fh_in)
                for row in reader:
                    filename = row["File"]
                    ref = row["Ref"]
                    output_prefix = row["OutputPrefix"]
                    if UmberjackWork.is_sam(filename):
                        yield self.sam_ref_window_iter(samfile=filename, ref=ref, output_prefix=output_prefix)
                    elif UmberjackWork.is_msa_fasta(filename):
                        yield self.msa_fasta_window_iter(msa_fasta=filename, ref=ref, output_prefix=output_prefix)


    def input_csv_iter(self):
        """
        Parses self.input_csv to get list of input samfiles/multiple sequence aligned fasta files
        and their corresponding output files.
        Does not do any error checking.  Use UmberjackWork.check_input_csv() for error checking.
        Must have headers 'File,Ref,OutputPrefix'.  Should be comma-delimited.
        :return iterator of UmberjackWork.Job namedtuples: For each (samfile, ref, output prefix) or (msa fasta, ref, output prefix) combination,
            creates a new Job namedtuple
        """
        if self.input_csv:
            with open(self.input_csv, 'rU') as fh_in:
                reader = csv.DictReader(fh_in)
                for row in reader:
                    filename = row["File"]
                    ref = row["Ref"]
                    output_prefix = row["OutputPrefix"]
                    if UmberjackWork.is_msa_fasta(filename):
                        yield UmberjackWork.Job(File=filename, Ref=ref, OutputPrefix=output_prefix, FileType=UmberjackWork.InputFileTypes.MSA)
                    else:
                        yield UmberjackWork.Job(File=filename, Ref=ref, OutputPrefix=output_prefix, FileType=UmberjackWork.InputFileTypes.SAM)



    @staticmethod
    def is_msa_fasta(filename):
        """
        Returns if file is a fasta.
        Since there is no standard file suffix for fastas, simply checks if there is a line starting with ">" in the file
        and a non-zero length sequence.
        :param str filename:  filepath to file
        :return bool:  whether the file is probably a fasta
        """
        if not os.path.exists(filename) or not os.path.getsize(filename):
            raise ValueError("File " + filename + " does not exist or is empty")

        first_header = Utility.get_first_header(filename)
        seqlen = Utility.get_len_1st_seq(filename)
        if first_header and seqlen:
            return True

        return  False


    @staticmethod
    def is_sam(filename):
        """
        Returns if the file is a sam.
        Does a quick check for headers and a tab-delimited line with at least 11 fields.
        :return:
        """
        if not os.path.exists(filename) or not os.path.getsize(filename):
            raise ValueError("File " + filename + " does not exist or is empty")

        has_header = sam.sam_handler.has_header(filename)
        has_samfields =  sam.sam_handler.has_samfields(filename)
        return has_header and has_samfields


    @staticmethod
    def check_sam(sam_filename, ref=None):
        """
        Checks that either sam_filename is defined and points to valid SAM file(s).
        Checks that all sam file has header and is queryname sorted.
        Checks that if ref isn't specified, then there are refs in the sam header.

        :param str sam_filename: filepath to sam
        :param str ref:  reference contig name in sam or None if all refs specified in sam header should be checked
        :return bool:  True if everything passes
        :raises ValueError:  if there are issues
        """
        if not os.path.exists(sam_filename) or not os.path.getsize(sam_filename):
            raise ValueError("SAM file " + sam_filename + " does not exist or is empty")

        if not sam_handler.is_query_sort(sam_filename):
            raise ValueError("SAM file " + sam_filename + " is not queryname sorted")

        refs = sam_handler.get_refs(sam_filename)
        if not refs:
            raise ValueError("SAM file " +  sam_filename + " does not have references in header")
        elif ref and ref not in refs:
            raise ValueError("SAM file " +  sam_filename + " does not contain " + ref + " in its references")

        return True


    @staticmethod
    def check_msa_fasta(msa_fasta):
        """
        Checks that msa_fasta exists and is multiple sequence aligned.
        Checks that all sam file has header and is queryname sorted.
        Checks that if ref isn't specified, then there are refs in the sam header.

        :return bool:  True if everything passes
        :raises ValueError:  if there are issues
        """
        if not os.path.exists(msa_fasta) or not os.path.getsize(msa_fasta):
            raise ValueError("Multiple Sequence Aligned fasta file " + msa_fasta + " does not exist or is empty")

        name = None
        seq = None
        last_seq_len = None
        names = set([])
        with open(msa_fasta, 'rU') as fh_in:
            for line in fh_in:
                if line.startswith(">"):
                    name = line[1:].split()[0]  # Get the read ID and ignore the description
                    formatted_name = Utility.newick_nice_name(name)
                    if formatted_name in names:
                        raise ValueError("Sequences with duplicate names after replacing periods, semicolons, colons, brackets with underscore")
                    names.add(formatted_name)

                    if seq is not None:
                        if last_seq_len is not None and last_seq_len != len(seq):
                            raise ValueError("Fasta is not multiple sequence aligned. Every sequence should be same length")
                        last_seq_len = len(seq)

                    seq = ""

                elif name is None:
                    raise ValueError("Comments and other non-sequence related lines are not allowed in multiple sequence aligned fasta")
                else:
                    seq += line.rstrip()

            # handle last line
            if seq == "":
                raise ValueError("Missing sequence for last entry")
            elif last_seq_len is not None and last_seq_len != len(seq):
                raise ValueError("Fasta is not multiple sequence aligned. Every sequence should be same length")

        return True


    def check_input(self):
        """
        Checks that one of sam_filename or input_csv is specified.
        Checks that all sam files specified have headers and are queryname sorted.
        Checks that if ref isn't specified, then ref is stored in each specified sam file.
        Checks that input CSV has correct headers 'File,Ref,OutputPrefix'.
        Checks that the files specified in the input csv are either sam files or multiple sequence aligned fasta files.

        :return bool:  True if everything passes
        :raises ValueError:  if there are issues
        """
        if not self.sam_filename and not self.msa_fasta and not self.input_csv:
            raise ValueError("Must specify one of sam_filename, msa_fasta, or input_csv")
        elif self.sam_filename:
            UmberjackWork.check_sam(self.sam_filename, self.ref)
        elif self.msa_fasta:
            UmberjackWork.check_msa_fasta(self.msa_fasta)
        elif self.input_csv:
            self.check_input_csv()

        return True


    @staticmethod
    def create_full_msa_fasta(sam_filename, output_prefix, ref, mapping_cutoff, read_qual_cutoff,
                          is_insert, is_mask_stop_codon):
        """
        Creates a pseudo multiple-sequence aligned fasta file for all reads using pairwise alignment from a SAM file.
        Does not filter based on breadth thresholds or N's.  But it does mask low quality bases and conflicts.

        :return: path to multiple sequence aligned fasta file of all reads
        :rtype : str
        :param str sam_filename: filepath to sam file of read alignments to the reference.  Must have header and must be queryname sorted.
        :param str output_prefix: output file prefix including parent directory but excluding file suffix
        :param str ref: reference name
        :param int mapping_cutoff: minimum mapping quality cutoff.  Reads aligned with map quality below this are thrown out.
        :param int read_qual_cutoff: read quality cutoff.  Bases below this cutoff are converted to N's.
        :param bool is_insert:  If True, then keeps insertions with respect to reference
        :param bool is_mask_stop_codon:  If True, then masks stop codons
        """

        msa_fasta_filename = output_prefix + ".msa.fasta"

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
    def create_dup_tsv(sam_filename, output_prefix, ref, mapping_cutoff, read_qual_cutoff, is_insert):
        """
        Spits the read duplicates to tab separated file

        :return: path to multiple sequence aligned fasta file of all reads
        :rtype : str
        :param str sam_filename: filepath to sam file of read alignments to the reference.  Must have header and must be queryname sorted.
        :param str output_prefix: output file prefix including parent directory but excluding file suffix
        :param str ref: reference name
        :param int mapping_cutoff: minimum mapping quality cutoff.  Reads aligned with map quality below this are thrown out.
        :param int read_qual_cutoff: read quality cutoff.  Bases below this cutoff are converted to N's.
        :param bool is_insert:  If True, then keeps insertions with respect to reference
        """

        dup_tsv_filename =output_prefix + ".dup.tsv"

        LOGGER.debug("Start Duplicate Reads TSV " + dup_tsv_filename + " from SAM for ref " + ref)
        if not os.path.exists(dup_tsv_filename) or os.path.getsize(dup_tsv_filename) <= 0:
            total_uniq = sam.sam_handler.write_dup_record_tsv(sam_filename=sam_filename, ref= ref,
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
        with open(listfile, 'rU') as fh_in:
            for line in fh_in:
                line = line.rstrip()
                items.extend([line])
        return items


    def sam_ref_iter(self):
        """
        Yields an iterator for each sample-ref combination
        :return iter(str, str): Iterator for (samfile, ref) (samfiletuple
        """
        queue_samfiles = []
        if self.sam_filename:
            queue_samfiles.extend([self.sam_filename])

        for queue_next_samfile in queue_samfiles:
            if self.ref:
                queue_refs = [self.ref]
            else:
                queue_refs = sam_handler.get_refs(queue_next_samfile)

            for queue_next_ref in queue_refs:
                yield queue_next_samfile, queue_next_ref


    @staticmethod
    def get_auto_outdir(pardir, input_file, ref):
        """
        Returns the auto-generated output directory path of <pardir>/<input_file file basename without file suffix>/<ref>
        This is the auto-generated output directory if none is given to Umberjack.

        :param str pardir: filepath to parent dir
        :param str input_file:  filepath to sam file or msa fasta file
        :param str ref:  name of reference.  If None, then outputs to <pardir>/<input_file file basename without file suffix>
        :return str: filepath in which Umberjack will spit out all the intermediary and end results
        """
        basename = os.path.basename(input_file)
        file_prefix = os.path.splitext(basename)[0]
        if ref:
            outdir = pardir + os.sep + file_prefix + os.sep + ref
        else:
            outdir = pardir + os.sep + file_prefix

        return outdir


    @staticmethod
    def get_auto_outprefix(pardir, input_file, ref):
        """
        Returns the auto-generated output filename prefix of
        <pardir>/<input_file file basename without file suffix>/<ref>/<input_file file basename without file suffix>.<ref>

        This is the auto-generated output filename prefix if none is given to Umberjack.

        :param str pardir: filepath to parent dir
        :param str input_file:  filepath to sam file or msa fasta file
        :param str ref:  name of reference.  If None, then outputs to <pardir>/<input_file file basename without file suffix>
        :return str: filepath in which Umberjack will spit out all the intermediary and end results
        """
        basename = os.path.basename(input_file)
        file_prefix = os.path.splitext(basename)[0]
        if ref:
            return pardir + os.sep + file_prefix + os.sep + ref + os.sep + file_prefix + "." + ref
        else:
            return pardir + os.sep + file_prefix + os.sep + file_prefix


    def sam_ref_window_iter(self, samfile, ref, output_prefix=None):
        """
        Yields an iterator over windows generated from sam & ref input.
        :param str samfile:  filepath to sam file.
        :param str ref:  reference config name.  If None, then will iterate over all references specified in sam header.
        :param str output_prefix:  filepath to output files, excluding the file suffix
        :returns iterator of dict: for every window in the sam-ref combination, yields a dict of window args that can be passed to eval_window()
        """
        if not output_prefix:
            output_prefix = UmberjackWork.get_auto_outprefix(pardir=self.out_dir, input_file=samfile, ref=ref)


        if self.debug:
            UmberjackWork.create_full_msa_fasta(sam_filename=samfile, output_prefix=output_prefix, ref=ref,
                                  mapping_cutoff=self.map_qual_cutoff, read_qual_cutoff=self.read_qual_cutoff,
                                  is_insert=self.insert, is_mask_stop_codon=self.mask_stop_codon)

        if self.debug and self.remove_duplicates:
            UmberjackWork.create_dup_tsv(sam_filename=samfile, output_prefix=output_prefix, ref=ref,
                           mapping_cutoff=self.map_qual_cutoff, read_qual_cutoff=self.read_qual_cutoff,
                           is_insert=self.insert)

        # All nucleotide positions are 1-based
        if not self.start_nucpos:
            queue_ref_start_nucpos = 1
        else:
            queue_ref_start_nucpos = self.start_nucpos
        if not self.end_nucpos:
            queue_ref_end_nucpos = sam.sam_handler.get_reflen(samfile, ref)
        else:
            queue_ref_end_nucpos =  self.end_nucpos


        total_windows = int(math.ceil(((queue_ref_end_nucpos - queue_ref_start_nucpos + 1) - self.window_size + 1)/self.window_slide))
        LOGGER.debug("Queuing " + str(total_windows) + " total windows for sam=" + samfile + " ref=" + ref)


        for start_window_nucpos in range(queue_ref_start_nucpos, queue_ref_end_nucpos-self.window_size+2, self.window_slide):
            end_window_nucpos = start_window_nucpos + self.window_size - 1

            window_args = {"window_depth_cutoff": self.window_depth_cutoff,
                           "window_breadth_cutoff": self.window_breadth_cutoff,
                           "start_window_nucpos": start_window_nucpos,
                           "end_window_nucpos": end_window_nucpos,
                           "ref": ref,
                           "output_prefix": output_prefix,
                           "sam_filename": samfile,
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


    def msa_fasta_window_iter(self, msa_fasta, ref=None, output_prefix=None):
        """
        :param str msa_fasta:  filepath to multiple sequence aligned (MSA) fasta file.
        :param str ref:  reference config name.  Only used to name output directory if output_prefix is not given.
        :param str output_prefix:  filepath to output files, excluding the file suffix.  If not given,
            then uses <out_dir>/<msa file basename prefix>/<ref> as the output directory.
        :returns iterator of dict: for every window in the sam-ref combination, yields a dict of window args that can be passed to eval_window()
        """

        if not self.start_nucpos:
            queue_ref_start_nucpos = 1
        else:
            queue_ref_start_nucpos = self.start_nucpos

        if not self.end_nucpos:
            queue_ref_end_nucpos = Utility.get_len_1st_seq(msa_fasta)
        else:
            queue_ref_end_nucpos =  self.end_nucpos

        if not output_prefix:
            output_prefix = UmberjackWork.get_auto_outprefix(pardir=self.out_dir, input_file=msa_fasta, ref=ref)

        # TODO:  remove duplicates from msa fasta

        # All nucleotide positions are 1-based
        total_windows = int(math.ceil(((queue_ref_end_nucpos - queue_ref_start_nucpos + 1) - self.window_size + 1)/self.window_slide))
        LOGGER.debug("Queuing " + str(total_windows) + " total windows for msa fasta=" + msa_fasta)
        LOGGER.debug("queue_ref_end_nucpos=" + str(queue_ref_end_nucpos))
        LOGGER.debug("queue_ref_end_nucpos-self.window_size+1=" + str(queue_ref_end_nucpos-self.window_size+1))

        for start_window_nucpos in range(queue_ref_start_nucpos, queue_ref_end_nucpos-self.window_size+2, self.window_slide):
            end_window_nucpos = start_window_nucpos + self.window_size - 1
            window_args = {"window_depth_cutoff": self.window_depth_cutoff,
                           "window_breadth_cutoff": self.window_breadth_cutoff,
                           "start_window_nucpos": start_window_nucpos,
                           "end_window_nucpos": end_window_nucpos,
                           "output_prefix": output_prefix,
                           "msa_fasta": msa_fasta,
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


    def window_iter(self):
        """
        Yields an iterator for each samfile-ref-window, or msa_fasta window.

        Creates the full MSA fasta and duplicates file for every samfile-ref.
        :return:
        """
        if self.sam_filename:
            return self.sam_ref_window_iter(samfile=self.sam_filename, ref=self.ref, output_prefix=None)
        elif self.msa_fasta:
            return self.msa_fasta_window_iter(msa_fasta=self.msa_fasta, ref=self.ref, output_prefix=None)
        elif self.input_csv:
            return self.input_csv_window_iter()
        else:
            raise ValueError("Nothing to iterate over.  Must define sam_filename or input_csv")


    def job_iter(self):
        """
        Iterates through all the (sam, ref, output prefix) or (msa fasta, ref, output prefix) combinations.
        :return iterator of Job namedtuples:
        """
        if self.sam_filename:
            for samfile, ref in self.sam_ref_iter():
                output_prefix = UmberjackWork.get_auto_outprefix(pardir=self.out_dir, input_file=samfile, ref=ref)
                yield UmberjackWork.Job(File=samfile, Ref=ref, OutputPrefix=output_prefix, FileType=UmberjackWork.InputFileTypes.SAM)

        elif self.msa_fasta:
            output_prefix = UmberjackWork.get_auto_outprefix(pardir=self.out_dir, input_file=self.msa_fasta, ref=self.ref)
            yield UmberjackWork.Job(File=self.msa_fasta, Ref=self.ref, OutputPrefix=output_prefix, FileType=UmberjackWork.InputFileTypes.MSA)

        elif self.input_csv:
            yield self.input_csv_iter()


    def tabulate_results(self):
        """
        Tabulates results for a samfile & ref into a single csv file.
        Plots results if applicable.
        """
        LOGGER.debug("About to tabulate results")
        if self.mode == MODE_DNDS:
            kwds_list = []
            for job in self.job_iter():
                output_csv_filename = job.OutputPrefix + ".dnds.csv"
                out_dir = os.path.dirname(output_csv_filename)

                if job.FileType == UmberjackWork.InputFileTypes.SAM:
                    ref_len = sam_handler.get_reflen(sam_filename=job.File, ref=job.Ref)
                else:
                    ref_len = Utility.get_len_1st_seq(job.File)

                kwds_list.append({"dnds_tsv_dir": out_dir,
                                  "ref": job.Ref,
                                  "ref_nuc_len": ref_len,
                                  "output_csv_filename": output_csv_filename
                })

            self.pool.spread_work(func=slice_miseq.tabulate_dnds, work_args_iter=kwds_list)

            plot_kws_list = []
            for job in self.job_iter():
                output_csv_filename = job.OutputPrefix + ".dnds.csv"
                plot_kws_list.append({"dnds_csv": output_csv_filename})

            self.pool.spread_work(func=plotter.plot_dnds, work_args_iter=plot_kws_list)


        elif self.mode == MODE_GTR_RATE:
            kwds_list = []
            for job in self.job_iter():
                output_csv_filename = job.OutputPrefix + ".gtr.csv"
                out_dir = os.path.dirname(output_csv_filename)

                kwds_list.append({"fasttree_output_dir": out_dir,
                                  "output_csv_filename": output_csv_filename
                })

            self.pool.spread_work(func=slice_miseq.tabulate_rates, work_args_iter=kwds_list)

        LOGGER.debug("Done tabulating results")


    def eval_windows(self):
        """
        spread_work window processes asynchronously.
        :return:
        """
        LOGGER.debug("About to spread window work")
        self.pool.spread_work(func=eval_window, work_args_iter=self.window_iter())
        LOGGER.debug("Done spread window work")


    def start(self):
        """
        Main entry method for the pipeline.
        :return:
        """

        self.setup()

        # spread_work windows for every sam/msa fasta, ref into the same pool
        self.eval_windows()


        # Tabulate results for every sam/msa fasta, ref
        self.tabulate_results()


    def stop(self):
        """
        Kill all processes
        :return:
        """
        self.pool.terminate()


