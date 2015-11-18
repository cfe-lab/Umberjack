#!/usr/bin/python
"""
Main entry point for Umberjack.
"""
import sys
import logging
import argparse

from config_arg_parse import ConfigArgParser
import config.settings
import UmberjackWork as UmberjackPool
import traceback

LOGGER = logging.getLogger(__name__)



def main():
    """
    Parses commandline arguments and kicks off mpi or multiprocessing versions of window evaluations.
    """
    umberjackworker = None
    try:
        parser = ConfigArgParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description="Version " + str(config.settings.VERSION),
                                 fromfile_prefix_chars='@')
        parser.add_argument("-f", help="Full filepath to config file.  If defined, then ignores all other commandline arguments.")
        parser.add_argument("--sam_filename", default=None,
                            help=("Full filepath to SAM alignment file.  Must be queryname sorted and must have header." +
                                  "One of --sam_filename, --msa_fasta, or --input_csv must be specified."))
        parser.add_argument("--msa_fasta", default=None,
                            help=("Full filepath to Multiple Sequence Aligned (MSA) fasta file." +
                                  "One of --sam_filename, --msa_fasta, or --input_csv must be specified."))
        parser.add_argument("--ref", default="",
                            help=("Name of reference contig.  " +
                                  "If --sam_filename is supplied, --ref must match one of the referenced specified in the SAM header.  " +
                                  "If --sam_filename is supplied but --ref is not supplied, "
                                  "then Umberjack processes every reference specified in the SAM file header."))
        parser.add_argument("--out_dir", default=config.settings.DEFAULT_OUT_DIR,
                            help="Output directory in which the pipeline will write all its results and intermediate files.  "
                                 "Output files will follow format " +
                                 "<out_dir>/<sam_filename or msa_fasta prefix>/<ref>/<sam_filename or msa_fasta prefix>.<ref>.*")
        parser.add_argument("--input_csv", default=None,
                            help="Filepath to comma-separated values (CSV) file containing table of " +
                                 "input sam files or input Multiple Sequence Aligned (MSA) fasta files and their corresponding output file prefixes.  " +
                                 "Headers should be 'File,Ref,OutputPrefix'. " +
                                 "'File' field should be the filepath of the input sam or MSA fasta file.  " +
                                 "'Ref' field should be the name of the reference contig in the sam file or that the MSA fasta file is aligned to.  " +
                                 "'OutputPrefix' field should include any parent directory and the output filename without any file suffix.  "
                                 "Checks if input file is a sam file or a multiple-sequence aligned fasta file from its file suffix and content.  " +
                                 "This option is to allow for processing of many sam or MSA fasta files at once.")
        parser.add_argument("--map_qual_cutoff", type=int, default=config.settings.DEFAULT_MAP_QUAL_CUTOFF,
                            help="Mapping quality threshold below which alignments are ignored.")
        parser.add_argument("--read_qual_cutoff", type=int, default=config.settings.DEFAULT_READ_QUAL_CUTOFF,
                            help="Read quality threshold below which bases are converted to Ns.")
        parser.add_argument("--max_prop_n", type=float, default=config.settings.DEFAULT_MAX_PROP_N,
                            help="Maximum fraction of Ns allowed in single-end reads or merged paired reads below which the read is ignored.")
        parser.add_argument("--window_size", type=int, default=config.settings.DEFAULT_WINDOW_SIZE,
                            help="Window size in nucleotides.")
        parser.add_argument("--window_slide", type=int, default=config.settings.DEFAULT_WINDOW_SLIDE,
                            help="Number of bases to slide each window by.")
        parser.add_argument("--window_breadth_cutoff", type=float, default=config.settings.DEFAULT_WINDOW_BREADTH_CUTOFF,
                            help="fraction of window that a single-end or merged paired read must cover with non-gap and non-N"
                                 " nucleotides.  Below this threshold, the read is omitted from the window.")
        parser.add_argument("--window_depth_cutoff", type=int, default=config.settings.DEFAULT_WINDOW_DEPTH_CUTOFF,
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
        parser.add_argument("--threads_per_window", type=int, default=config.settings.DEFAULT_THREADS_PER_WINDOW,
                            help="threads allotted per window.")
        parser.add_argument("--concurrent_windows", type=int, default=config.settings.DEFAULT_CONCURRENT_WINDOWS,
                            help="Max number of windows to process concurrently. Ignored when --mpi is defined.")
        parser.add_argument("--hyphy_exe", default=config.settings.DEFAULT_HYPHY_EXE,
                            help="full filepath of HYPHYMP executable.  Default: taken from PATH")
        parser.add_argument("--hyphy_basedir", default=config.settings.DEFAULT_HYPHY_BASEDIR,
                            help="full filepath of custom HyPhy template batch files.")
        parser.add_argument("--hyphy_libdir", default=config.settings.DEFAULT_HYPHY_LIBDIR,
                            help="full filepath of HyPhy base directory of out-of-box template batch files.")
        parser.add_argument("--fastree_exe", default=config.settings.DEFAULT_FASTTREEMP_EXE,
                            help="full filepath of FastTreeMP or FastTree executable.  Default: taken from PATH")
        parser.add_argument("--mode", default=config.settings.DEFAULT_MODE, choices=[UmberjackPool.MODE_DNDS, UmberjackPool.MODE_GTR_RATE, UmberjackPool.MODE_COUNT_SUBS],
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
                LOGGER.info("Using config file " + args.f + ".  Ignoring all other commandline arguments.")
            args = parser.parse_args([parser.fromfile_prefix_chars +  args.f])

        check_single = [x for x in [args.sam_filename, args.msa_fasta,  args.input_csv] if x]
        if len(check_single) != 1:
            errmsg = "Must specify one of --sam_filename, --msa_fasta, or --input_csv arguments"
            print errmsg
            LOGGER.error(errmsg)
            parser.print_help()
            sys.exit()

        # deep copy of arguments excluding empty values
        args_msg = ""
        eval_windows_args = {}
        for key, val in vars(args).iteritems():
            args_msg += "{}={} ".format(key, val)
            if val is not None and key != "f":
                eval_windows_args[key] = val
        LOGGER.info("Launched umberjack.py with settings: " + args_msg)

        # Automatically set mask_stop_codon to True for dN/dS analysis.
        if args.mode == UmberjackPool.MODE_DNDS:
            eval_windows_args["mask_stop_codon"] = True
            eval_windows_args["remove_duplicates"] = True
            LOGGER.warn("Auto-setting --mask_stop_codon to True since Umberjack is in " + UmberjackPool.MODE_DNDS + " mode")
            LOGGER.warn("Auto-setting --remove_duplicates to True since Umberjack is in " + UmberjackPool.MODE_DNDS + " mode")

        # If the user has mpi4py installed and user specified -mpi in commandline,
        # then runs the MPI version otherwise runs the multiprocessing version on current node
        if args.mpi:
            try:
                import pool.MPIPool
                # Ignore the concurrent_windows commandline arg and uses the number of processors indicated by mpirun command
                eval_windows_args.pop("concurrent_windows", None)
                LOGGER.info("Running MPI Version. Ignoring --concurrent_windows flag.  Using mpirun node arguments.")
                processpool = pool.MPIPool.MPIPool()
            except ImportError:
                LOGGER.warn("You must install mpi4py module in order to leverage multiple nodes.  Running on single node.")
                eval_windows_args.pop("mpi", None)
                import pool.OneNodePool
                processpool = pool.OneNodePool.OneNodePool(args.concurrent_windows)
        else:
            import pool.OneNodePool
            processpool = pool.OneNodePool.OneNodePool(args.concurrent_windows)

        eval_windows_args["pool"] = processpool


        umberjackworker = UmberjackPool.UmberjackWork(**eval_windows_args)
        umberjackworker.start()

    except:
        LOGGER.exception("Umberjack failure.")
        sys.exit(1)
    finally:
        if umberjackworker:
            umberjackworker.stop()  # Kill all processes


if __name__ == '__main__':
    config.settings.setup_logging()
    main()
