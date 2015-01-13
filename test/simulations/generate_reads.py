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

PICARD_BIN_DIR = sys.argv[7]

BOWTIE_OUT_DIR = sys.argv[8]
PROCS = int(sys.argv[9])
SEED = int(sys.argv[10])


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
     "-l",  "250", # length of read
     #"-f", "2", # fold coverage
     "-f", "2", # fold coverage
     "-m",  "346", # mean fragment size
     "-s",  "75", # std dev fragment size
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
if not os.path.exists(BOWTIE_OUT_DIR):
    os.makedirs(BOWTIE_OUT_DIR)


bowtie_output_prefix = BOWTIE_OUT_DIR + os.sep + os.path.basename(art_output_prefix)
bowtie_log = bowtie_output_prefix + ".bowtie.log"
bowtie_db_prefix = BOWTIE_OUT_DIR + os.sep + os.path.basename(consensus_fasta).replace(".fasta", "")
BOWTIE_BUILD_CMD = ["bowtie2-build",
                    consensus_fasta,  # fasta to align against
                    bowtie_db_prefix]  # prefix of db index

# Super high gap penalty, very sensitive local alignment
BOWTIE_CMD = ["bowtie2",
              "--local",
              "-D", "30",  # max allowed failed seed extensions
              "-R", "4", # max reseeds for repetitive seeds
              "-N", "1", # max mismatches during seed alignment
              "-L", "10", # seed length
              "-i", "S,1,0.25",  # interval between seed strings.  Interval(readlen) = 1 + 0.25 * sqrt(readlen) = 4.9
              "-x", bowtie_db_prefix,  # reference index
              "-t",  # output timing
              "-1", art_output_prefix + ".1.fq",
              "-2", art_output_prefix + ".2.fq",
              "-S", bowtie_output_prefix + ".consensus.bowtie.sam",  # samfile output
              "--rdg", "100,3",  # read affine gap open, gap extension penalty
              "--rfg", "100,3",  # ref affine gap open, gap extension penalty
              "--phred33",   # sanger phred scoring
              "-p", str(PROCS)]

# Super high gap penalty, very sensitive local alignment
#--very-sensitive-local
#Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
# -D 30 -R 4 -N 1 -L 10 -i S,1,0.25
BOWTIE_ERRFREE_CMD = ["bowtie2",
              "--local",
              "-D", "30",  # max allowed failed seed extensions
              "-R", "4", # max reseeds for repetitive seeds
              "-N", "1", # max mismatches during seed alignment
              "-L", "10", # seed length
              "-i", "S,1,0.25",  # interval between seed strings.  Interval(readlen) = 1 + 0.25 * sqrt(readlen)
              "-x", bowtie_db_prefix,  # reference index
              "-t",  # output timing
              "-1", art_output_prefix + ".errFree.1.fq",
              "-2", art_output_prefix + ".errFree.2.fq",
              "-S", bowtie_output_prefix + ".errFree.consensus.bowtie.sam",  # samfile output
              "--rdg", "100,3",  # read affine gap open, gap extension penalty
              "--rfg", "100,3",  # ref affine gap open, gap extension penalty
              "--phred33",   # sanger phred scoring
              "-p", str(PROCS)]

with open(bowtie_log, 'w') as fh_log:
    print "Logging to " + bowtie_log
    print "About to execute " + " ".join(BOWTIE_BUILD_CMD)
    subprocess.check_call(BOWTIE_BUILD_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

    print "About to execute " + " ".join(BOWTIE_CMD)
    subprocess.check_call(BOWTIE_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)

    print "About to execute " + " ".join(BOWTIE_ERRFREE_CMD)
    subprocess.check_call(BOWTIE_ERRFREE_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)




# Get the coverage stats and histo
picard_logfile = bowtie_output_prefix + ".picard.wgsmetrics.log"
print "Logging to " + picard_logfile
if os.path.exists(picard_logfile):
    os.remove(picard_logfile)
for input_sam, aln_ref_fasta in [(art_output_prefix + ".sam", reference_fasta),
                                 (art_output_prefix + ".errFree.sam", reference_fasta),
                                 (bowtie_output_prefix + ".consensus.bowtie.sam", consensus_fasta),
                                 (bowtie_output_prefix + ".errFree.consensus.bowtie.sam", consensus_fasta)]:

    output_sam = input_sam.replace(".sam", ".sort.sam")
    output_wgsmetrics = input_sam.replace(".sam", ".picard.wgsmetrics")
    PICARD_SORT_SAM_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "picard.jar",
                  "SortSam",
                  "INPUT=" + input_sam,
                  "OUTPUT=" + output_sam,
                  "SORT_ORDER=coordinate"]
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



# Get per-bp coverage along reference consensus.  This table is passed to R for plotting.
samtools_logfile = bowtie_output_prefix + ".samtools.log"
print "Logging to " + samtools_logfile
if os.path.exists(samtools_logfile):
    os.remove(samtools_logfile)
for input_sam in [art_output_prefix + ".sort.sam",
                  art_output_prefix + ".errFree.sort.sam",
                  bowtie_output_prefix + ".consensus.bowtie.sort.sam",
                  bowtie_output_prefix + ".errFree.consensus.bowtie.sort.sam"]:

    output_bam = input_sam.replace(".sam", ".bam")
    SAMTOOLS_BAM_CMD = ["samtools", "view",
                        "-S",  # input is Sam format
                        "-b", # output to bam format
                        "-o", output_bam, # output bam file,
                        input_sam ] # input samfile
    output_depth_tsv = input_sam.replace(".sam", ".cov.tsv")
    SAMTOOLS_DEPTH_CMD = ["samtools", "depth",
                          "-q", "20",  # min base quality 20
                          "-Q", "20",  # min mapping quality 20
                          output_bam ]
    with open(samtools_logfile, 'a') as fh_log,  open(output_depth_tsv, 'w') as fh_out_cov:
        print "About to execute " + " ".join(SAMTOOLS_BAM_CMD)
        subprocess.check_call(SAMTOOLS_BAM_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)
        print "About to execute " + " ".join(SAMTOOLS_DEPTH_CMD)
        subprocess.check_call(SAMTOOLS_DEPTH_CMD, env=os.environ, stdout=fh_out_cov, stderr=fh_log)