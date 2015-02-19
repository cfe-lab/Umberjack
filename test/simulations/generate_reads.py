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
                outrow["Conserve"] = aln.get_conserve(nuc_pos_0based, is_count_ambig=True, is_count_gaps=True, is_count_pad=False)
                outrow["Entropy"] = aln.get_shannon_entropy(nuc_pos_0based, is_count_ambig=True, is_count_gaps=True, is_count_pad=False)
                outrow["NucDepth"] = aln.get_depth(nuc_pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                codon_pos_0based = nuc_pos_0based / Utility.NUC_PER_CODON
                outrow["CodonDepth"] = codon_depths[codon_pos_0based]
                writer.writerow(outrow)


# Get stats on merged read error, length, N's
for msa_fasta in [reference_fasta,      # full population, no sequencing
                      art_output_prefix + ".msa.fasta",
                      art_output_prefix + ".errFree.msa.fasta",
                      bwa_output_prefix + ".consensus.bwa.msa.fasta",
                      bwa_output_prefix + ".errFree.consensus.bwa.msa.fasta"]:
        nuc_conserve_csv = msa_fasta.replace(".fasta", ".conserve.csv").replace(".msa", "")


import Bio.SeqIO as SeqIO
import re
full_popn_recdict = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))
consec_rec = SeqIO.read(consensus_fasta, "fasta")

output_cmp_msa_csv  = art_output_prefix + ".cmp.msa.csv"
output_err_msa_csv  = art_output_prefix + ".err.msa.csv"
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

    for source, msa_fasta in [("Orig", art_output_prefix + ".msa.fasta"),
                              ("OrigErrFree", art_output_prefix + ".errFree.msa.fasta"),
                              ("Aln", bwa_output_prefix + ".consensus.bwa.msa.fasta"),
                              ("AlnErrFree", bwa_output_prefix + ".errFree.consensus.bwa.msa.fasta")]:

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

Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R")
# Output the R config file (workaround cuz can't seem to set commandline parameters for rmd knitr scripts)
Rcov_config_file = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R" + os.sep + "small_cov.config")
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
    fh_out_rconfig.write("CMP_READ_ERR_CSV={}".format(output_cmp_msa_csv) + "\n")

Rscript_cmd = ("library(knitr); " +
               "setwd('{}'); ".format(Rscript_wdir) +
               "spin('small_cov.R', knit=FALSE); " +
               "knit2html('small_cov.Rmd', stylesheet='markdown_bigwidth.css')")
subprocess.check_call(["Rscript", "-e", Rscript_cmd], env=os.environ)
shutil.copy(Rscript_wdir + os.sep + "small_cov.html", bwa_output_prefix + ".cov.html")







