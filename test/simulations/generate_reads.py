# ASSUMES that art_illumina is the PATH env var
import os
import subprocess
import shutil
import Utility
import csv
import sam.sam_handler as sam_handler
import sys
import logging

# If fragments are 500bp and read pairs do not overlap:  for 1x coverage, 10k indiv * 9kbp per indiv / 500bp per pair = 180000 pairs
# If fragments are 400bp and read pairs do not overlap: for 1x coverage, 10k indiv * 9kbp per indiv / 400bp per pair = 225000 pairs.
# If fragments are 500bp and read pairs overlap by 100bp:  225000 pairs x 500bp per pair / (10kindiv * 9kbp) = 1.25x overall population coverage.



LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)



ART_BIN_DIR = sys.argv[1]
ART_QUAL_PROFILE_TSV1 = sys.argv[2]
ART_QUAL_PROFILE_TSV2 = sys.argv[3]
reference_fasta = sys.argv[4]
consensus_fasta = sys.argv[5]
art_output_prefix = sys.argv[6]
ART_FOLD_COVER = sys.argv[7]
art_mean_insert_size = sys.argv[8]
art_stddev_insert_size = sys.argv[9]

PICARD_BIN_DIR = sys.argv[10]

BWA_OUT_DIR = sys.argv[11]
PROCS = int(sys.argv[12])
SEED = int(sys.argv[13])
INDELIBLE_RATES_CSV = sys.argv[14]

CONSENSUS_NAME = "consensus"
MIN_BASE_Q = 20
MIN_MAP_Q = 20
MAX_PROP_N = 0.1

# Simulate reads with ART
#############################################

if os.path.exists(art_output_prefix + ".sam"):
    LOGGER.warn("Not regenerating ART simulated reads")
else:

    if not os.path.exists(os.path.dirname(art_output_prefix)):
        os.makedirs(os.path.dirname(art_output_prefix))

    ART_CMD=[ART_BIN_DIR + os.sep + "art_illumina",
         "-na", # don't output alignment file
         "-ef", # create both error-free and normal reads
         "-sam", # create sam alignment
         "-p",  # paired end,
         "-rs", str(SEED),
         "-1", ART_QUAL_PROFILE_TSV1, # 1st read quality  profile
         "-2", ART_QUAL_PROFILE_TSV2,  # 2nd read quality profile
         "-d", "read", # read id prefix
         "-i",  reference_fasta, # dna reference fasta
         "-ir", "0", # 1st read insertion rate
         "-ir2",  "0", # 2nd read insertion rate
         "-dr",  "0", # 1st read deletion rate
         "-dr2",  "0", # 2nd read deletion rate
         "-qs", "2",  # Bump up the quality scores of every base in 1st mate so that average error rate = 0.006
         "-qs2", "2",  # Bump up the quality scores of every base in 2nd mate so that average error rate = 0.006
         "-l",  "250", # length of read
         "-f", ART_FOLD_COVER, # fold coverage
         "-m",  art_mean_insert_size, # mean fragment size
         "-s",  art_stddev_insert_size, # std dev fragment size
         "-o",  art_output_prefix  # output prefix
         ]

    logfile = art_output_prefix + ".art.log"
    with open(logfile, 'w') as fh_log:
        print "Logging to " + logfile
        print "About to execute " + " ".join(ART_CMD)
        subprocess.check_call(ART_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

    # ART creates fastq files with suffix 1.fq.  We want .1.fq suffix.  Rename files.
    shutil.move(art_output_prefix + "1.fq", art_output_prefix + ".1.fq")
    shutil.move(art_output_prefix + "2.fq", art_output_prefix + ".2.fq")
    # ART creates error free sam with suffix "_errFree.sam".  We want .errFree.sam.  Rename files
    shutil.move(art_output_prefix + "_errFree.sam", art_output_prefix + ".errFree.sam")

    # ART creates an error free sam.  Generate error free fq from error free sam
    PICARD_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                  "SamToFastq",
                  "INPUT=" + art_output_prefix + ".errFree.sam",
                  "FASTQ=" + art_output_prefix + ".errFree.1.fq",
                  "SECOND_END_FASTQ=" + art_output_prefix + ".errFree.2.fq"]
    picard_logfile = art_output_prefix + ".picard.sam2fastq.log"
    with open(picard_logfile, 'w') as fh_log:
        print "Logging to " + picard_logfile
        print "About to execute " + " ".join(PICARD_CMD)
        subprocess.check_call(PICARD_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)



# Align reads against population consensus
#############################################

bwa_output_prefix = BWA_OUT_DIR + os.sep + os.path.basename(art_output_prefix)
if os.path.exists(bwa_output_prefix + ".consensus.bwa.sam"):
    LOGGER.warn("Not realigning reads")
else:
    if not os.path.exists(BWA_OUT_DIR):
        os.makedirs(BWA_OUT_DIR)

    bwa_log = bwa_output_prefix + ".bwa.log"
    bwa_db_prefix = BWA_OUT_DIR + os.sep + os.path.basename(consensus_fasta).replace(".fasta", "")
    BWA_BUILD_CMD = ["/home/thuy/programs/bwa/bwa-0.7.8/bwa", "index",
                     "-p", bwa_db_prefix, # prefix of db index
                        consensus_fasta,  # fasta to align against
                        ]

    # bwa seems to perform better than bowtie for aligning the very divergent sequences
    with open(bwa_log, 'w') as fh_log:
        print "Logging to " + bwa_log
        print "About to execute " + " ".join(BWA_BUILD_CMD)
        subprocess.check_call(BWA_BUILD_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

        for fq_prefix, output_sam in [(art_output_prefix, bwa_output_prefix + ".consensus.bwa.sam"),
                               (art_output_prefix + ".errFree", bwa_output_prefix + ".errFree.consensus.bwa.sam")]:
            BWA_CMD = ["/home/thuy/programs/bwa/bwa-0.7.8/bwa",
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
                print "About to execute " + " ".join(BWA_CMD) + " output to " + output_sam
                subprocess.check_call(BWA_CMD, env=os.environ, stdout=fh_out_sam, stderr=fh_log)



# Get Coverage Stats
####################################################################

if os.path.exists(art_output_prefix + ".cov.tsv"):
    LOGGER.warn("Not regenerating coverage stats")
else:

    # Get the coverage stats and histo
    # Get alignment stats
    picard_logfile = bwa_output_prefix + ".picard.metrics.log"
    print "Logging to " + picard_logfile
    if os.path.exists(picard_logfile):
        os.remove(picard_logfile)
    for input_sam, aln_ref_fasta in [(art_output_prefix + ".sam", reference_fasta),
                                     (art_output_prefix + ".errFree.sam", reference_fasta),
                                     (bwa_output_prefix + ".consensus.bwa.sam", consensus_fasta),
                                     (bwa_output_prefix + ".errFree.consensus.bwa.sam", consensus_fasta)]:

        output_sam = input_sam.replace(".sam", ".sort.sam")
        output_query_sam = input_sam.replace(".sam", ".sort.query.sam")
        output_wgsmetrics = input_sam.replace(".sam", ".picard.wgsmetrics")
        PICARD_SORT_SAM_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                      "SortSam",
                      "INPUT=" + input_sam,
                      "OUTPUT=" + output_sam,
                      "SORT_ORDER=coordinate"]
        # Sort by queryname for sliding_window_tree.py
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
                      "MAX_INSERT_SIZE=" + str(int(art_mean_insert_size) + 2 * int(art_stddev_insert_size)),   # max insert size after which pair is considered  chimeric
                      "METRIC_ACCUMULATION_LEVEL=ALL_READS" # get metrics for all reads
                      ]
        with open(picard_logfile, 'a') as fh_log:
            print "About to execute " + " ".join(PICARD_SORT_SAM_CMD)
            subprocess.check_call(PICARD_SORT_SAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

            print "About to execute " + " ".join(PICARD_WGS_METRICS_CMD)
            subprocess.check_call(PICARD_WGS_METRICS_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

            print "About to execute " + " ".join(PICARD_SORT_QUERYNAME_SAM_CMD)
            subprocess.check_call(PICARD_SORT_QUERYNAME_SAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

            print "About to execute " + " ".join(PICARD_ALN_METRICS_CMD)
            subprocess.check_call(PICARD_ALN_METRICS_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)




    # Get per-bp coverage along reference consensus.  This table is passed to R for plotting.
    samtools_logfile = bwa_output_prefix + ".samtools.log"
    print "Logging to " + samtools_logfile
    if os.path.exists(samtools_logfile):
        os.remove(samtools_logfile)
    for input_sam in [art_output_prefix + ".sort.sam",
                      art_output_prefix + ".errFree.sort.sam",
                      bwa_output_prefix + ".consensus.bwa.sort.sam",
                      bwa_output_prefix + ".errFree.consensus.bwa.sort.sam"]:

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
            print "About to execute " + " ".join(SAMTOOLS_BAM_CMD)
            subprocess.check_call(SAMTOOLS_BAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)
            print "About to execute " + " ".join(SAMTOOLS_DEPTH_CMD) + " to " + output_depth_tsv
            subprocess.check_call(SAMTOOLS_DEPTH_CMD, env=os.environ, stdout=fh_out_cov, stderr=fh_log)

# Get Conservation And Entropy Stats for Multiple Sequence Alignment
####################################################################

if os.path.exists(art_output_prefix + ".conserve.csv"):
    LOGGER.warn("Not regenerating conservation & entropy stats for MSA")
else:
    # Create a multiple sequence alignment for the full population and entire genome from the aligned reads.
    # We will use these for obtaining conservation and entropy at each position for the aligned reads.
    consensus_len = Utility.get_seq2len(consensus_fasta)[CONSENSUS_NAME]
    for sam, ref in [(art_output_prefix + ".sort.query.sam", None),
                     (art_output_prefix + ".errFree.sort.query.sam", None),
                     (bwa_output_prefix + ".consensus.bwa.sort.query.sam", CONSENSUS_NAME),
                     (bwa_output_prefix + ".errFree.consensus.bwa.sort.query.sam", CONSENSUS_NAME)]:
        msa_fasta = sam.replace(".sort.query.sam", ".msa.fasta")
        sam_handler.create_msa_slice_from_sam(sam_filename=sam,
                                              ref=ref,
                                              out_fasta_filename=msa_fasta,
                                              mapping_cutoff=MIN_MAP_Q,
                                              read_qual_cutoff=MIN_BASE_Q,
                                              max_prop_N=MAX_PROP_N,
                                              breadth_thresh=0.0,
                                              start_pos=None, end_pos=None,  # slice is the entire length of the genome
                                              is_insert=False, ref_len=consensus_len)



    # Get conservation and entropy stats at each nucleotide position for the full population, full simulated reads, aligned reads
    for msa_fasta in [reference_fasta,      # full population, no sequencing
                      art_output_prefix + ".msa.fasta",
                      art_output_prefix + ".errFree.msa.fasta",
                      bwa_output_prefix + ".consensus.bwa.msa.fasta",
                      bwa_output_prefix + ".errFree.consensus.bwa.msa.fasta"]:
        nuc_conserve_csv = msa_fasta.replace(".fasta", ".conserve.csv").replace(".msa", "")

        with open(nuc_conserve_csv, 'w') as fh_out:
            LOGGER.debug("About to output conservation, entropy from " + msa_fasta + " to " + nuc_conserve_csv)
            writer = csv.DictWriter(fh_out, fieldnames=["NucSite", "Conserve", "Entropy", "NucDepth", "CodonDepth"])
            writer.writeheader()

            aln = Utility.Consensus()
            aln.parse(msa_fasta)
            codon_depths = Utility.get_total_codons_by_pos(msa_fasta)  # TODO:  cleaner if we put this in Consensus class

            for nuc_pos_0based in range(0, aln.get_alignment_len()):
                outrow = dict()
                outrow["NucSite"] = nuc_pos_0based + 1
                outrow["Conserve"] = aln.get_conserve(nuc_pos_0based)
                outrow["Entropy"] = aln.get_metric_entropy(nuc_pos_0based)
                outrow["NucDepth"] = aln.get_depth(nuc_pos_0based)
                codon_pos_0based = nuc_pos_0based / Utility.NUC_PER_CODON
                outrow["CodonDepth"] = codon_depths[codon_pos_0based]
                writer.writerow(outrow)




# Pass Info to R for plotting into HTML format
####################################################################

Rscript_wdir =  os.path.abspath(os.path.dirname(__file__) + os.sep + "R")
# Output the R config file (workaround cuz can't seem to set commandline parameters for rmd knitr scripts)
Rcov_config_file = os.path.abspath(os.path.dirname(__file__) + os.sep + "R" + os.sep + "small_cov.config")
with open(Rcov_config_file, 'w') as fh_out_rconfig:
    fh_out_rconfig.write("FULL_POPN_CONSERVE_CSV={}".format(reference_fasta.replace(".fasta", ".conserve.csv")) + "\n")
    fh_out_rconfig.write("ART_FOLD_COVER={}".format(ART_FOLD_COVER) + "\n")
    fh_out_rconfig.write("ORIG_ERR_FREE_COV_TSV={}".format(art_output_prefix + ".errFree.cov.tsv") + "\n")
    fh_out_rconfig.write("ALN_ERR_FREE_COV_TSV={}".format(bwa_output_prefix + ".errFree.consensus.bwa.cov.tsv") + "\n")
    fh_out_rconfig.write("ORIG_ERR_FREE_WGS_METRICS={}".format(art_output_prefix + ".errFree.picard.wgsmetrics") + "\n")
    fh_out_rconfig.write("ALN_ERR_FREE_WGS_METRICS={}".format(bwa_output_prefix + ".errFree.consensus.bwa.picard.wgsmetrics") + "\n")
    fh_out_rconfig.write("ORIG_ERR_FREE_CONSERVE_CSV={}".format(art_output_prefix + ".errFree.conserve.csv") + "\n")
    fh_out_rconfig.write("ALN_ERR_FREE_CONSERVE_CSV={}".format(bwa_output_prefix + ".errFree.consensus.bwa.conserve.csv") + "\n")
    fh_out_rconfig.write("ORIG_COV_TSV={}".format(art_output_prefix + ".cov.tsv") + "\n")
    fh_out_rconfig.write("ALN_COV_TSV={}".format(bwa_output_prefix + ".consensus.bwa.cov.tsv") + "\n")
    fh_out_rconfig.write("ORIG_WGS_METRICS={}".format(art_output_prefix + ".picard.wgsmetrics") + "\n")
    fh_out_rconfig.write("ALN_WGS_METRICS={}".format(bwa_output_prefix + ".consensus.bwa.picard.wgsmetrics") + "\n")
    fh_out_rconfig.write("ORIG_CONSERVE_CSV={}".format(art_output_prefix + ".conserve.csv") + "\n")
    fh_out_rconfig.write("ALN_CONSERVE_CSV={}".format(bwa_output_prefix + ".consensus.bwa.conserve.csv") + "\n")
    fh_out_rconfig.write("INDELIBLE_RATES_CSV={}".format(INDELIBLE_RATES_CSV) + "\n")

Rscript_cmd = ("library(knitr); " +
               "setwd('{}'); ".format(Rscript_wdir) +
               "spin('small_cov.R', knit=FALSE); " +
               "knit2html('small_cov.Rmd', stylesheet='markdown_bigwidth.css')")
subprocess.check_call(["Rscript", "-e", Rscript_cmd], env=os.environ)
shutil.copy(Rscript_wdir + os.sep + "small_cov.html", bwa_output_prefix + ".cov.html")


# At each nuc site, count number of errors in Original Reads
# At each nuc site, count number of conserved base -> non-conserved base due to sequencing error in Original Reads
# At each nuc site, count number of mutation -> conserved-base due to sequencing error in Original Reads
# At each nuc site, count number of errors in Aligned Reads
# At each nuc site, count number of conserved base -> non-conserved base due to sequencing error in Aligned Reads
# At each nuc site, count number of mutation -> conserved-base due to sequencing error in Aligned Reads
