# ASSUMES that art_illumina is the PATH env var
import os
import subprocess
import shutil


# If fragments are 500bp and read pairs do not overlap:  for 1x coverage, 10k indiv * 9kbp per indiv / 500bp per pair = 180000 pairs
# If fragments are 400bp and read pairs do not overlap: for 1x coverage, 10k indiv * 9kbp per indiv / 400bp per pair = 225000 pairs.
# If fragments are 500bp and read pairs overlap by 100bp:  225000 pairs x 500bp per pair / (10kindiv * 9kbp) = 1.25x overall population coverage.

import sys

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
NUM_INDIV = sys.argv[14]
NUM_CODON_SITES = sys.argv[15]


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
if not os.path.exists(BWA_OUT_DIR):
    os.makedirs(BWA_OUT_DIR)


bwa_output_prefix = BWA_OUT_DIR + os.sep + os.path.basename(art_output_prefix)
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
            print "About to execute " + " ".join(BWA_CMD)
            subprocess.check_call(BWA_CMD, env=os.environ, stdout=fh_out_sam, stderr=fh_log)



# Get the coverage stats and histo
picard_logfile = bwa_output_prefix + ".picard.wgsmetrics.log"
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

    with open(picard_logfile, 'a') as fh_log:
        print "About to execute " + " ".join(PICARD_SORT_SAM_CMD)
        subprocess.check_call(PICARD_SORT_SAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)
        print "About to execute " + " ".join(PICARD_WGS_METRICS_CMD)
        subprocess.check_call(PICARD_WGS_METRICS_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)
        print "About to execute " + " ".join(PICARD_SORT_QUERYNAME_SAM_CMD)
        subprocess.check_call(PICARD_SORT_QUERYNAME_SAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)



# Get per-bp coverage along reference consensus.  This table is passed to R for plotting.
samtools_logfile = bwa_output_prefix + ".samtools.log"
MIN_BASE_Q = 20
MIN_MAP_Q = 5
print "Logging to " + samtools_logfile
if os.path.exists(samtools_logfile):
    os.remove(samtools_logfile)
for input_sam in [art_output_prefix + ".sort.sam",
                  art_output_prefix + ".errFree.sort.sam",
                  bwa_output_prefix + ".consensus.bwa.sort.sam",
                  bwa_output_prefix + ".errFree.consensus.bwa.sort.sam"]:

    output_bam = input_sam.replace(".sam", ".bam")
    SAMTOOLS_BAM_CMD = ["samtools", "view",
                        "-S",  # input is Sam format
                        "-b", # output to bam format
                        "-o", output_bam, # output bam file,
                        input_sam ] # input samfile
    output_depth_tsv = input_sam.replace(".sam", ".cov.tsv")
    SAMTOOLS_DEPTH_CMD = ["samtools", "depth",
                          "-q", str(MIN_BASE_Q),  # min base quality
                          "-Q", str(MIN_MAP_Q),  # min mapping quality 20
                          output_bam ]
    with open(samtools_logfile, 'a') as fh_log,  open(output_depth_tsv, 'w') as fh_out_cov:
        print "About to execute " + " ".join(SAMTOOLS_BAM_CMD)
        subprocess.check_call(SAMTOOLS_BAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)
        print "About to execute " + " ".join(SAMTOOLS_DEPTH_CMD)
        subprocess.check_call(SAMTOOLS_DEPTH_CMD, env=os.environ, stdout=fh_out_cov, stderr=fh_log)

# Plot the coverage and qualties into HTML
Rscript_wdir =  os.path.abspath(os.path.dirname(__file__) + os.sep + "R")
# Output the R config file (workaround cuz can't seem to set commandline parameters for rmd knitr scripts)
Rcov_config_file = os.path.abspath(os.path.dirname(__file__) + os.sep + "R" + os.sep + "small_cov.config")
with open(Rcov_config_file, 'w') as fh_out_rconfig:
    fh_out_rconfig.write("NUM_INDIV={}".format(NUM_INDIV) + "\n")
    fh_out_rconfig.write("NUM_CODON_SITES={}".format(NUM_CODON_SITES) + "\n")
    fh_out_rconfig.write("ART_FOLD_COVER={}".format(ART_FOLD_COVER) + "\n")
    fh_out_rconfig.write("ORIG_ERR_FREE_COV_TSV={}".format(art_output_prefix + ".errFree.sort.cov.tsv",) + "\n")
    fh_out_rconfig.write("ALN_ERR_FREE_COV_TSV={}".format(bwa_output_prefix + ".errFree.consensus.bwa.sort.cov.tsv",) + "\n")
    fh_out_rconfig.write("ORIG_ERR_FREE_WGS_METRICS={}".format(art_output_prefix + ".errFree.picard.wgsmetrics",) + "\n")
    fh_out_rconfig.write("ALN_ERR_FREE_WGS_METRICS={}".format(bwa_output_prefix + ".errFree.consensus.bwa.picard.wgsmetrics",) + "\n")
    fh_out_rconfig.write("ORIG_COV_TSV={}".format(art_output_prefix + ".sort.cov.tsv",) + "\n")
    fh_out_rconfig.write("ALN_COV_TSV={}".format(bwa_output_prefix + ".consensus.bwa.sort.cov.tsv",) + "\n")
    fh_out_rconfig.write("ORIG_WGS_METRICS={}".format(art_output_prefix + ".picard.wgsmetrics",) + "\n")
    fh_out_rconfig.write("ALN_WGS_METRICS={}".format(bwa_output_prefix + ".consensus.bwa.picard.wgsmetrics",) + "\n")
Rscript_cmd = ("library(knitr); " +
               "setwd('{}'); ".format(Rscript_wdir) +
               "spin('small_cov.R')")
subprocess.check_call(["Rscript", "-e", Rscript_cmd], env=os.environ)
shutil.copy(Rscript_wdir + os.sep + "small_cov.html", bwa_output_prefix + ".cov.html")
