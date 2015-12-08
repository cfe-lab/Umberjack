"""
Simulates paired-end read sequencing with ART simulator.
Aligns the reads to the reference with BWA.
Converts sam to reads with Picard.
Outputs coverage stats to HTML using R Knitr package.
"""
import os
import subprocess
import shutil
import Utility
import csv
import sam.sam_handler as sam_handler
import sys
import logging
import config.settings as settings
import ConfigParser
import sim_helper
# If fragments are 500bp and read pairs do not overlap:  for 1x coverage, 10k indiv * 9kbp per indiv / 500bp per pair = 180000 pairs
# If fragments are 400bp and read pairs do not overlap: for 1x coverage, 10k indiv * 9kbp per indiv / 400bp per pair = 225000 pairs.
# If fragments are 500bp and read pairs overlap by 100bp:  225000 pairs x 500bp per pair / (10kindiv * 9kbp) = 1.25x overall population coverage.


settings.setup_logging()
LOGGER = logging.getLogger(__name__)


SECTION = "sim"

config_file = sys.argv[1]
config = ConfigParser.RawConfigParser()
config.read(config_file)

OUTDIR = os.path.dirname(config_file)  # Output directory for simulated data

SEED = config.getint(SECTION, "SEED")
FILENAME_PREFIX = config.get(SECTION, "FILENAME_PREFIX")


ART_BIN_DIR = sim_helper.get_path_str(config.get(SECTION, "ART_BIN_DIR"), OUTDIR)
ART_QUAL_PROFILE_TSV1 = sim_helper.get_path_str(config.get(SECTION, "ART_QUAL_PROFILE_TSV1"), OUTDIR)
ART_QUAL_PROFILE_TSV2 = sim_helper.get_path_str(config.get(SECTION, "ART_QUAL_PROFILE_TSV2"), OUTDIR)
ART_FOLD_COVER = config.getfloat(SECTION, "ART_FOLD_COVER")
ART_MEAN_FRAG = config.getint(SECTION, "ART_MEAN_FRAG")
ART_STDEV_FRAG = config.getint(SECTION, "ART_STDEV_FRAG")
ART_INSERT_RATE1 = config.getfloat(SECTION, "ART_INSERT_RATE1")
ART_INSERT_RATE2 = config.getfloat(SECTION, "ART_INSERT_RATE2")
ART_DEL_RATE1 = config.getfloat(SECTION, "ART_DEL_RATE1")
ART_DEL_RATE2 = config.getfloat(SECTION, "ART_DEL_RATE2")
ART_QUAL_SHIFT1 = config.getint(SECTION, "ART_QUAL_SHIFT1")
ART_QUAL_SHIFT2 = config.getint(SECTION, "ART_QUAL_SHIFT2")
ART_READ_LENGTH = config.getint(SECTION, "ART_READ_LENGTH")


PICARD_BIN_DIR = sim_helper.get_path_str(config.get(SECTION, "PICARD_BIN_DIR"), OUTDIR)
BWA_BIN_DIR = sim_helper.get_path_str(config.get(SECTION, "BWA_BIN_DIR"), OUTDIR)

PROCS = config.getint(SECTION, "PROCS")

# TODO:  don't hard code these - get them from config file
ART_OUTPUT_PREFIX = OUTDIR + os.sep +  "reads" + os.sep + FILENAME_PREFIX + ".reads"
REFERENCE_FASTA = OUTDIR + os.sep + "fullpopn" + os.sep +  FILENAME_PREFIX + "_TRUE.fasta"  # Use the INDELible genome fasta as the population reference
CONSENSUS_FASTA = OUTDIR + os.sep + "fullpopn" + os.sep +  FILENAME_PREFIX + ".consensus.fasta"
BWA_OUT_DIR = OUTDIR + os.sep +  "aln"
INDELIBLE_RATES_CSV = OUTDIR + os.sep + "fullpopn" + os.sep + FILENAME_PREFIX + "_RATES.csv"


CONSENSUS_NAME = "consensus"
MIN_BASE_Q = 20       # TODO:  get this from the config file
MIN_MAP_Q = 20
MAX_PROP_N = 0.1

# Simulate reads with ART
#############################################

if (os.path.exists(ART_OUTPUT_PREFIX + ".sam") and os.path.getsize(ART_OUTPUT_PREFIX + ".sam") and
        os.path.exists(ART_OUTPUT_PREFIX + ".reads.1.fq") and os.path.getsize(ART_OUTPUT_PREFIX + ".reads.1.fq") and
        os.path.exists(ART_OUTPUT_PREFIX + ".reads.2.fq") and os.path.getsize(ART_OUTPUT_PREFIX + ".reads.2.fq")):
    LOGGER.warn("Not regenerating ART simulated reads")
else:

    if not os.path.exists(os.path.dirname(ART_OUTPUT_PREFIX)):
        os.makedirs(os.path.dirname(ART_OUTPUT_PREFIX))

    ART_CMD=[ART_BIN_DIR + os.sep + "art_illumina",
         "-na", # don't output alignment file
         "-ef", # create both error-free and normal reads
         "-sam", # create sam alignment
         "-p",  # paired end,
         "-rs", str(SEED),
         "-1", ART_QUAL_PROFILE_TSV1, # 1st read quality  profile
         "-2", ART_QUAL_PROFILE_TSV2,  # 2nd read quality profile
         "-d", "read", # read id prefix
         "-i",  REFERENCE_FASTA, # dna reference fasta
         "-ir", str(ART_INSERT_RATE1), # 1st read insertion rate
         "-ir2",  str(ART_INSERT_RATE2), # 2nd read insertion rate
         "-dr",  str(ART_DEL_RATE1), # 1st read deletion rate
         "-dr2",  str(ART_DEL_RATE2), # 2nd read deletion rate
         "-qs", str(ART_QUAL_SHIFT1),  # Bump up the quality scores of every base in 1st mate so that average error rate = 0.006
         "-qs2", str(ART_QUAL_SHIFT2),  # Bump up the quality scores of every base in 2nd mate so that average error rate = 0.006
         "-l",  str(ART_READ_LENGTH), # length of read
         "-f", str(ART_FOLD_COVER), # fold coverage
         "-m",  str(ART_MEAN_FRAG), # mean fragment size
         "-s",  str(ART_STDEV_FRAG), # std dev fragment size
         "-o",  ART_OUTPUT_PREFIX  # output prefix
         ]

    logfile = ART_OUTPUT_PREFIX + ".art.log"
    with open(logfile, 'w') as fh_log:
        LOGGER.debug( "Logging to " + logfile)
        LOGGER.debug( "About to execute " + " ".join(ART_CMD))
        subprocess.check_call(ART_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

    # ART creates fastq files with suffix 1.fq.  We want .1.fq suffix.  Rename files.
    shutil.move(ART_OUTPUT_PREFIX + "1.fq", ART_OUTPUT_PREFIX + ".1.fq")
    shutil.move(ART_OUTPUT_PREFIX + "2.fq", ART_OUTPUT_PREFIX + ".2.fq")
    # ART creates error free sam with suffix "_errFree.sam".  We want .errFree.sam.  Rename files
    shutil.move(ART_OUTPUT_PREFIX + "_errFree.sam", ART_OUTPUT_PREFIX + ".errFree.sam")

    # ART creates an error free sam.  Generate error free fq from error free sam
    PICARD_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                  "SamToFastq",
                  "INPUT=" + ART_OUTPUT_PREFIX + ".errFree.sam",
                  "FASTQ=" + ART_OUTPUT_PREFIX + ".errFree.1.fq",
                  "SECOND_END_FASTQ=" + ART_OUTPUT_PREFIX + ".errFree.2.fq"]
    picard_logfile = ART_OUTPUT_PREFIX + ".picard.sam2fastq.log"
    with open(picard_logfile, 'w') as fh_log:
        LOGGER.debug( "Logging to " + picard_logfile)
        LOGGER.debug( "About to execute " + " ".join(PICARD_CMD))
        subprocess.check_call(PICARD_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)



# Align reads against population consensus
#############################################

bwa_output_prefix = BWA_OUT_DIR + os.sep + os.path.basename(ART_OUTPUT_PREFIX)
if os.path.exists(bwa_output_prefix + ".bwa.sam") and os.path.getsize(bwa_output_prefix + ".bwa.sam"):
    LOGGER.warn("Not realigning reads")
else:
    if not os.path.exists(BWA_OUT_DIR):
        os.makedirs(BWA_OUT_DIR)

    bwa_log = bwa_output_prefix + ".bwa.log"
    bwa_db_prefix = BWA_OUT_DIR + os.sep + os.path.basename(CONSENSUS_FASTA).replace(".fasta", "")
    BWA_BUILD_CMD = [BWA_BIN_DIR + os.sep + "bwa", "index",
                     "-p", bwa_db_prefix, # prefix of db index
                        CONSENSUS_FASTA,  # fasta to align against
                        ]

    # bwa seems to perform better than bowtie for aligning the very divergent sequences
    with open(bwa_log, 'w') as fh_log:
        LOGGER.debug( "Logging to " + bwa_log)
        LOGGER.debug( "About to execute " + " ".join(BWA_BUILD_CMD))
        subprocess.check_call(BWA_BUILD_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

        for fq_prefix, output_sam in [(ART_OUTPUT_PREFIX, bwa_output_prefix + ".bwa.sam"),
                               (ART_OUTPUT_PREFIX + ".errFree", bwa_output_prefix + ".errFree.bwa.sam")]:
            BWA_CMD = [BWA_BIN_DIR + os.sep + "bwa",
                  "mem",
                  "-t", "10", # threads
                  "-k", "5",  # seed len
                  "-d", "4800", # max extension
                  "-r", "1", # reseeding
                  "-B", "2", # mismatch penalty
                  "-T", "10", # min align score
                  bwa_db_prefix,  # reference index
                  fq_prefix + ".1.fq",
                  fq_prefix + ".2.fq"]
            with open(output_sam, 'w') as fh_out_sam:
                LOGGER.debug( "About to execute " + " ".join(BWA_CMD) + " output to " + output_sam)
                subprocess.check_call(BWA_CMD, env=os.environ, stdout=fh_out_sam, stderr=fh_log)



# Get Coverage Stats
####################################################################

if os.path.exists(ART_OUTPUT_PREFIX + ".cov.tsv") and os.path.getsize(ART_OUTPUT_PREFIX + ".cov.tsv"):
    LOGGER.warn("Not regenerating coverage stats")
else:

    # Get the coverage stats and histo
    # Get alignment stats
    picard_logfile = bwa_output_prefix + ".picard.metrics.log"
    LOGGER.debug( "Logging to " + picard_logfile)
    if os.path.exists(picard_logfile):
        os.remove(picard_logfile)
    for input_sam, aln_ref_fasta in [(ART_OUTPUT_PREFIX + ".sam", REFERENCE_FASTA),
                                     (ART_OUTPUT_PREFIX + ".errFree.sam", REFERENCE_FASTA),
                                     (bwa_output_prefix + ".bwa.sam", CONSENSUS_FASTA),
                                     (bwa_output_prefix + ".errFree.bwa.sam", CONSENSUS_FASTA)]:

        output_sam = input_sam.replace(".sam", ".sort.sam")
        output_query_sam = input_sam.replace(".sam", ".sort.query.sam")
        output_wgsmetrics = input_sam.replace(".sam", ".picard.wgsmetrics")
        PICARD_SORT_SAM_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                      "SortSam",
                      "INPUT=" + input_sam,
                      "OUTPUT=" + output_sam,
                      "SORT_ORDER=coordinate"]
        # Sort by queryname for umberjack.py
        PICARD_SORT_QUERYNAME_SAM_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                      "SortSam",
                      "INPUT=" + input_sam,
                      "OUTPUT=" + output_query_sam,
                      "SORT_ORDER=queryname"]

        # CollectWgsMetrics Needs a sam/bam sorted by coordinates
        PICARD_WGS_METRICS_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                      "CollectWgsMetrics",
                      "INPUT=" + output_sam,
                      "OUTPUT=" + output_wgsmetrics,
                      "REFERENCE_SEQUENCE=" + aln_ref_fasta,
                      "INCLUDE_BQ_HISTOGRAM=true"]

        # CollectAlignmentSummaryMetrics
        output_alnmetrics = input_sam.replace(".sam", ".picard.alnmetrics")
        PICARD_ALN_METRICS_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                      "CollectAlignmentSummaryMetrics",
                      "INPUT=" + output_sam,
                      "OUTPUT=" + output_alnmetrics,
                      "ASSUME_SORTED=false",  # read the sam header to determine if it's been sorted or not
                      "ADAPTER_SEQUENCE=null", # the ART generated reads don't have adapters
                      "REFERENCE_SEQUENCE=" + aln_ref_fasta,  # reference fasta.  If this is empty, then it doesn't collect metrics
                      "MAX_INSERT_SIZE=" + str(int(ART_MEAN_FRAG) + 2 * int(ART_STDEV_FRAG)),   # max insert size after which pair is considered  chimeric
                      "METRIC_ACCUMULATION_LEVEL=ALL_READS" # get metrics for all reads
                      ]
        with open(picard_logfile, 'a') as fh_log:
            LOGGER.debug( "About to execute " + " ".join(PICARD_SORT_SAM_CMD))
            subprocess.check_call(PICARD_SORT_SAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

            LOGGER.debug(  "About to execute " + " ".join(PICARD_WGS_METRICS_CMD))
            subprocess.check_call(PICARD_WGS_METRICS_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

            LOGGER.debug(  "About to execute " + " ".join(PICARD_SORT_QUERYNAME_SAM_CMD))
            subprocess.check_call(PICARD_SORT_QUERYNAME_SAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

            LOGGER.debug(  "About to execute " + " ".join(PICARD_ALN_METRICS_CMD))
            subprocess.check_call(PICARD_ALN_METRICS_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)




    # Get per-bp coverage along reference consensus.  This table is passed to R for plotting.
    samtools_logfile = bwa_output_prefix + ".samtools.log"
    LOGGER.debug(  "Logging to " + samtools_logfile)
    if os.path.exists(samtools_logfile):
        os.remove(samtools_logfile)
    for input_sam in [ART_OUTPUT_PREFIX + ".sort.sam",
                      ART_OUTPUT_PREFIX + ".errFree.sort.sam",
                      bwa_output_prefix + ".bwa.sort.sam",
                      bwa_output_prefix + ".errFree.bwa.sort.sam"]:

        output_bam = input_sam.replace(".sort.sam", ".sort.bam")
        SAMTOOLS_BAM_CMD = ["samtools", "view",
                            "-S",  # input is Sam format
                            "-b", # output to bam format
                            "-o", output_bam, # output bam file,
                            input_sam ] # input samfile
        output_depth_tsv = input_sam.replace(".sort.sam", ".cov.tsv")
        SAMTOOLS_DEPTH_CMD = ["samtools", "depth",
                              "-q", str(MIN_BASE_Q),  # min base quality
                              "-Q", str(MIN_MAP_Q),  # min mapping quality 20
                              output_bam ]
        with open(samtools_logfile, 'a') as fh_log,  open(output_depth_tsv, 'w') as fh_out_cov:
            LOGGER.debug( "About to execute " + " ".join(SAMTOOLS_BAM_CMD))
            subprocess.check_call(SAMTOOLS_BAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)
            LOGGER.debug( "About to execute " + " ".join(SAMTOOLS_DEPTH_CMD) + " to " + output_depth_tsv)
            subprocess.check_call(SAMTOOLS_DEPTH_CMD, env=os.environ, stdout=fh_out_cov, stderr=fh_log)

# Get Conservation And Entropy Stats for Multiple Sequence Alignment
####################################################################

if os.path.exists(ART_OUTPUT_PREFIX + ".conserve.csv") and os.path.getsize(ART_OUTPUT_PREFIX + ".conserve.csv"):
    LOGGER.warn("Not regenerating conservation & entropy stats for MSA")
else:
    # Create a multiple sequence alignment for the full population and entire genome from the aligned reads.
    # We will use these for obtaining conservation and entropy at each position for the aligned reads.
    consensus_len = Utility.get_seq2len(CONSENSUS_FASTA)[CONSENSUS_NAME]
    for sam, ref in [(ART_OUTPUT_PREFIX + ".sort.query.sam", None),
                     (ART_OUTPUT_PREFIX + ".errFree.sort.query.sam", None),
                     (bwa_output_prefix + ".bwa.sort.query.sam", CONSENSUS_NAME),
                     (bwa_output_prefix + ".errFree.bwa.sort.query.sam", CONSENSUS_NAME)]:
        msa_fasta = sam.replace(".sort.query.sam", ".msa.fasta")
        sam_handler.create_msa_slice_from_sam(sam_filename=sam, ref=ref, out_fasta_filename=msa_fasta,
                                              mapping_cutoff=MIN_MAP_Q, read_qual_cutoff=MIN_BASE_Q,
                                              max_prop_N=MAX_PROP_N, breadth_thresh=0.0, start_pos=0, end_pos=0,
                                              do_insert_wrt_ref=False, ref_len=consensus_len)



    # Get conservation and entropy stats at each nucleotide position for the full population, full simulated reads, aligned reads
    for msa_fasta in [REFERENCE_FASTA,      # full population, no sequencing
                      ART_OUTPUT_PREFIX + ".msa.fasta",
                      ART_OUTPUT_PREFIX + ".errFree.msa.fasta",
                      bwa_output_prefix + ".bwa.msa.fasta",
                      bwa_output_prefix + ".errFree.bwa.msa.fasta"]:
        nuc_conserve_csv = msa_fasta.replace(".fasta", ".conserve.csv").replace(".msa", "")

        with open(nuc_conserve_csv, 'w') as fh_out:
            LOGGER.debug("About to output conservation, entropy from " + msa_fasta + " to " + nuc_conserve_csv)
            writer = csv.DictWriter(fh_out, fieldnames=["CodonSite", "ConserveCodon", "EntropyCodon", "CodonDepth"])
            writer.writeheader()

            aln = Utility.Consensus()
            aln.parse(msa_fasta)

            for codon_pos_0based in range(0, aln.get_alignment_len()/3):
                outrow = dict()
                outrow["CodonSite"] = codon_pos_0based + 1
                outrow["ConserveCodon"] = aln.get_codon_conserve(codon_pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                outrow["EntropyCodon"] = aln.get_codon_shannon_entropy(codon_pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                outrow["CodonDepth"] = aln.get_codon_depth(codon_pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                writer.writerow(outrow)

output_cmp_msa_csv  = ART_OUTPUT_PREFIX + ".cmp.msa.csv"
output_err_msa_csv  = ART_OUTPUT_PREFIX + ".err.msa.csv"
if os.path.exists(output_cmp_msa_csv) and os.path.getsize(output_cmp_msa_csv) and os.path.exists(output_err_msa_csv) and os.path.getsize(output_err_msa_csv):
    LOGGER.warn("Not regenerating error csvs {} and {}".format(output_cmp_msa_csv, output_err_msa_csv))
else:
    # Get stats on merged read error, length, N's
    for msa_fasta in [REFERENCE_FASTA,      # full population, no sequencing
                          ART_OUTPUT_PREFIX + ".msa.fasta",
                          ART_OUTPUT_PREFIX + ".errFree.msa.fasta",
                          bwa_output_prefix + ".bwa.msa.fasta",
                          bwa_output_prefix + ".errFree.bwa.msa.fasta"]:
            nuc_conserve_csv = msa_fasta.replace(".fasta", ".conserve.csv").replace(".msa", "")


    import Bio.SeqIO as SeqIO
    import re
    full_popn_recdict = SeqIO.to_dict(SeqIO.parse(REFERENCE_FASTA, "fasta"))
    consec_rec = SeqIO.read(CONSENSUS_FASTA, "fasta")


    with open(output_cmp_msa_csv, 'w') as fh_out_cmp_msa_csv, open(output_err_msa_csv, 'w') as fh_out_err_msa_csv:
        writer_cmp = csv.DictWriter(fh_out_cmp_msa_csv,
                                fieldnames=["Read",     # Read name
                                            "Template",  # Read Template from Full Population
                                            "Source",  # Source of read  [Orig, OrigErrFree, Aln, AlnErrFree]
                                            "Start",  # 1-based start bp wrt ref
                                            "End",  # 1-basee end bp wrt ref
                                            "Ns",   # total N's
                                            "Consensus_Mismatch",  # Total read bases that don't match the consensus.  Tells us actual diversity.
                                            "Template_Consensus_Mismatch",  # Total bases in the template that don't match the full population consensus.  Tells us expected diversity.
                                            "SeqErr",   # Total read bases that don't match template.  Does not include N's
                                            "Conserve2Nonconserve", # Total template bases matching consensus converted to nonconsensus base in read.  Tells us if expected diversity is lower than actual
                                            "Nonconserve2Conserve",  # Total template bases that didn't match consensus convered to consensus base in read.  tells us if expected diversity is higher than actual.

                                ])
        writer_cmp.writeheader()

        writer_err = csv.DictWriter(fh_out_err_msa_csv,
                                fieldnames=["Read",     # Read name
                                            "Template",  # Read Template from Full Population
                                            "Source",  # Source of read  [Orig, OrigErrFree, Aln, AlnErrFree]
                                            "Start",  # 1-based start bp wrt ref
                                            "End",  # 1-basee end bp wrt ref
                                            "NucSite",
                                            "MutationType"  # N, Conserve2Nonconserve or Nonconserve2Conserve

                                ])
        writer_err.writeheader()

        for source, msa_fasta in [("Orig", ART_OUTPUT_PREFIX + ".msa.fasta"),
                                  ("OrigErrFree", ART_OUTPUT_PREFIX + ".errFree.msa.fasta"),
                                  ("Aln", bwa_output_prefix + ".bwa.msa.fasta"),
                                  ("AlnErrFree", bwa_output_prefix + ".errFree.bwa.msa.fasta")]:

            # Count the number of mismatches.
            # Count number of N's
            cmp_outrow = dict()
            cmp_outrow["Source"] = source

            for record in SeqIO.parse(msa_fasta, "fasta"):
                template, read = record.id.split("_")
                cmp_outrow["Template"] = template
                cmp_outrow["Read"] = read
                start = re.search(r"[^\-]", str(record.seq)).start() + 1  # 1-based
                end = re.search(r"[^\-][\-]*$", str(record.seq)).start() + 1  # 1-based
                cmp_outrow["Start"] = start
                cmp_outrow["End"] = end
                cmp_outrow["Ns"] = str(record.seq).count("N")

                # Only look non-left/right padded slice
                read_seq = str(record.seq[start-1:end])
                template_seq = str(full_popn_recdict[template].seq[start-1:end])
                consensus_seq = str(consec_rec.seq[start-1:end])

                ns = 0
                seq_err = 0
                consensus_mismatch = 0
                template_consensus_mismatch = 0
                conserve2Nonconserve = 0
                nonconserve2Conserve = 0
                err_outrow = dict()
                err_outrow["Source"] = source
                err_outrow["Template"] = template
                err_outrow["Read"] = read
                err_outrow["Start"] = start
                err_outrow["End"] = end
                for i in range(0, len(read_seq)):
                    if read_seq[i] != consensus_seq[i]:
                        consensus_mismatch += 1
                    if template_seq[i] != consensus_seq[i]:
                        template_consensus_mismatch += 1

                    if read_seq[i] != template_seq[i]:
                        err_outrow["NucSite"] = start + i
                        if read_seq[i] == "N":
                            ns += 1
                            err_outrow["MutationType"] = "N"
                        else:
                            seq_err += 1
                            if template_seq[i] == consensus_seq[i]:
                                conserve2Nonconserve += 1
                                err_outrow["MutationType"] = "Conserve2Nonconserve"
                            else:
                                nonconserve2Conserve += 1
                                err_outrow["MutationType"] = "Nonconserve2Conserve"

                        writer_err.writerow(err_outrow)

                cmp_outrow["Ns"] = ns
                cmp_outrow["SeqErr"] = seq_err
                cmp_outrow["Consensus_Mismatch"] = consensus_mismatch
                cmp_outrow["Conserve2Nonconserve"] = conserve2Nonconserve
                cmp_outrow["Nonconserve2Conserve"] = nonconserve2Conserve
                cmp_outrow["Template_Consensus_Mismatch"] = template_consensus_mismatch
                writer_cmp.writerow(cmp_outrow)




# Pass Info to R for plotting into HTML format
####################################################################
# if os.path.exists(bwa_output_prefix + ".cov.html") and os.path.getsize(bwa_output_prefix + ".cov.html"):
#     LOGGER.warn("Not regerating simulated data stats knitr html " + bwa_output_prefix + ".cov.html")
# else:
#     Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R")
#     # Output the R config file (workaround cuz can't seem to set commandline parameters for rmd knitr scripts)
#     Rcov_config_file = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R" + os.sep + "small_cov.config")
#     with open(Rcov_config_file, 'w') as fh_out_rconfig:
#         fh_out_rconfig.write("FULL_POPN_CONSERVE_CSV={}".format(REFERENCE_FASTA.replace(".fasta", ".conserve.csv")) + "\n")
#         fh_out_rconfig.write("ART_FOLD_COVER={}".format(ART_FOLD_COVER) + "\n")
#         fh_out_rconfig.write("ORIG_ERR_FREE_COV_TSV={}".format(ART_OUTPUT_PREFIX + ".errFree.cov.tsv") + "\n")
#         fh_out_rconfig.write("ALN_ERR_FREE_COV_TSV={}".format(bwa_output_prefix + ".errFree.bwa.cov.tsv") + "\n")
#         fh_out_rconfig.write("ORIG_ERR_FREE_WGS_METRICS={}".format(ART_OUTPUT_PREFIX + ".errFree.picard.wgsmetrics") + "\n")
#         fh_out_rconfig.write("ALN_ERR_FREE_WGS_METRICS={}".format(bwa_output_prefix + ".errFree.bwa.picard.wgsmetrics") + "\n")
#         fh_out_rconfig.write("ORIG_ERR_FREE_CONSERVE_CSV={}".format(ART_OUTPUT_PREFIX + ".errFree.conserve.csv") + "\n")
#         fh_out_rconfig.write("ALN_ERR_FREE_CONSERVE_CSV={}".format(bwa_output_prefix + ".errFree.bwa.conserve.csv") + "\n")
#         fh_out_rconfig.write("ORIG_COV_TSV={}".format(ART_OUTPUT_PREFIX + ".cov.tsv") + "\n")
#         fh_out_rconfig.write("ALN_COV_TSV={}".format(bwa_output_prefix + ".bwa.cov.tsv") + "\n")
#         fh_out_rconfig.write("ORIG_WGS_METRICS={}".format(ART_OUTPUT_PREFIX + ".picard.wgsmetrics") + "\n")
#         fh_out_rconfig.write("ALN_WGS_METRICS={}".format(bwa_output_prefix + ".bwa.picard.wgsmetrics") + "\n")
#         fh_out_rconfig.write("ORIG_CONSERVE_CSV={}".format(ART_OUTPUT_PREFIX + ".conserve.csv") + "\n")
#         fh_out_rconfig.write("ALN_CONSERVE_CSV={}".format(bwa_output_prefix + ".bwa.conserve.csv") + "\n")
#         fh_out_rconfig.write("INDELIBLE_RATES_CSV={}".format(INDELIBLE_RATES_CSV) + "\n")
#         fh_out_rconfig.write("CMP_READ_ERR_CSV={}".format(output_cmp_msa_csv) + "\n")
#
#     Rscript_cmd = ("library(knitr); " +
#                    "setwd('{}'); ".format(Rscript_wdir) +
#                    "spin('small_cov.R', knit=FALSE); " +
#                    "knit2html('small_cov.Rmd', stylesheet='markdown_bigwidth.css')")
#     subprocess.check_call(["Rscript", "-e", Rscript_cmd], env=os.environ)
#     shutil.copy(Rscript_wdir + os.sep + "small_cov.html", bwa_output_prefix + ".cov.html")







