# ASSUMES that art_illumina is the PATH env var
import os
import subprocess
import shutil


# If fragments are 500bp and read pairs do not overlap:  for 1x coverage, 10k indiv * 9kbp per indiv / 500bp per pair = 180000 pairs
# If fragments are 400bp and read pairs do not overlap: for 1x coverage, 10k indiv * 9kbp per indiv / 400bp per pair = 225000 pairs.
# If fragments are 500bp and read pairs overlap by 100bp:  225000 pairs x 500bp per pair / (10kindiv * 9kbp) = 1.25x overall population coverage.
#/home/thuy/programs/art/art_bin_VanillaIceCream/art_illumina -na -ef -sam -p -1 /home/thuy/programs/art/art_bin_VanillaIceCream/Illumina_profiles/EmpMiSeq250R1.txt -2 /home/thuy/programs/art/art_bin_VanillaIceCream/Illumina_profiles/EmpMiSeq250R2.txt -d read  -i /home/thuy/gitrepo/SlidingWindow/test/simulations/data/indelible/sample_genomes.100.fasta -ir 0 -ir2 0 -dr 0 -dr2 0 -l 250 -f 2 -m 400 -s 50 -o /home/thuy/gitrepo/SlidingWindow/test/simulations/data/indelible/cov1x/sample_genomes.100.1x &> /home/thuy/gitrepo/SlidingWindow/test/simulations/data/indelible/cov1x/sample_genomes.100.1x.art.log

# art_illumina [options] -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -o <outfile_prefix>
# art_illumina [options] -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -m <mean_fragsize> -s <std_fragsize> -o <outfile_prefix>
#
# ===== PARAMETERS =====
#
#   -1   --qprof1   the first-read quality profile
#   -2   --qprof2   the second-read quality profile
#   -amp --amplicon amplicon sequencing simulation
#   -d   --id       the prefix identification tag for read ID
#   -ef  --errfree  indicate to generate the zero sequencing errors SAM file as well the regular one
#                   NOTE: the reads in the zero-error SAM file have the same alignment positions
#                   as those in the regular SAM file, but have no sequencing errors
#   -f   --fcov     the fold of read coverage to be simulated or number of reads/read pairs generated for each amplicon
#   -h   --help     print out usage information
#   -i   --in       the filename of input DNA/RNA reference
#   -ir  --insRate  the first-read insertion rate (default: 0.00009)
#   -ir2 --insRate2 the second-read insertion rate (default: 0.00015)
#   -dr  --delRate  the first-read deletion rate (default:  0.00011)
#   -dr2 --delRate2 the second-read deletion rate (default: 0.00023)
#   -l   --len      the length of reads to be simulated
#   -m   --mflen    the mean size of DNA/RNA fragments for paired-end simulations
#   -mp  --matepair indicate a mate-pair read simulation
#   -nf  --maskN    the cutoff frequency of 'N' in a window size of the read length for masking genomic regions
#                   NOTE: default: '-nf 1' to mask all regions with 'N'. Use '-nf 0' to turn off masking
#   -na  --noALN    do not output ALN alignment file
#   -o   --out      the prefix of output filename
#   -p   --paired   indicate a paired-end read simulation or to generate reads from both ends of amplicons
#   -q   --quiet    turn off end of run summary
#   -qs  --qShift   the amount to shift every first-read quality score by
#   -qs2 --qShift2  the amount to shift every second-read quality score by
#                   NOTE: For -qs/-qs2 option, a positive number will shift up quality scores (the max is 93)
#                   that reduce substitution sequencing errors and a negative number will shift down
#                   quality scores that increase sequencing errors. If shifting scores by x, the error
#                   rate will be 1/(10^(x/10)) of the default profile.
#   -rs  --rndSeed  the seed for random number generator (default: system time in second)
#                   NOTE: using a fixed seed to generate two identical datasets from different runs
#   -s   --sdev     the standard deviation of DNA/RNA fragment size for paired-end simulations.
#   -sam --samout   indicate to generate SAM alignment file
#   -sp  --sepProf  indicate to use separate quality profiles for different bases (ATGC)
#                   NOTE: art will automatically switch to a mate-pair simulation if the given mean fragment size >= 2000
#
# ===== NOTES =====
#
# * ART by default selects a built-in quality score profile according to the read length specified for the run.
#
# * For single-end simulation, ART requires input sequence file, outputfile prefix, read length and fold coverage.
#
# * For paired-end simulation (except for amplicon sequencing), ART also requires the parameter values of
#   the mean and standard deviation of DNA/RNA fragment lengths

ART_BIN_DIR = "/home/thuy/programs/art/art_bin_VanillaIceCream"
DEFAULT_QUAL_PROFILE_TSV1 = "/home/thuy/programs/art/art_bin_VanillaIceCream/Illumina_profiles/EmpMiSeq250R1.txt"
DEFAULT_QUAL_PROFILE_TSV2 = "/home/thuy/programs/art/art_bin_VanillaIceCream/Illumina_profiles/EmpMiSeq250R2.txt"
#quality_profile_tsv = os.path.abspath("../data/quality_profiles")
reference_fasta = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small/small_population.mixedmut.fasta"
consensus_fasta = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small/small_population.mixedmut.consensus.fasta"
output_prefix = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small/small_population.mixedmut.reads"

PICARD_BIN_DIR = "/home/thuy/programs/picard/picard-tools-1.19"
PROCS = 10

ART_CMD=[ART_BIN_DIR + os.sep + "art_illumina",
     "-na", # don't output alignment file
     "-ef", # create both error-free and normal reads
     "-sam", # create sam alignment
     "-p",  # paired end
     "-1", DEFAULT_QUAL_PROFILE_TSV1, # 1st read quality  profile
     "-2", DEFAULT_QUAL_PROFILE_TSV2,  # 2nd read quality profile
     "-d", "read", # read id prefix
     "-i",  reference_fasta, # dna reference fasta
     "-ir", "0", # 1st read insertion rate
     "-ir2",  "0", # 2nd read insertion rate
     "-dr",  "0", # 1st read deletion rate
     "-dr2",  "0", # 2nd read deletion rate
     "-l",  "250", # length of read
     "-f", "2", # fold coverage
     "-m",  "346", # mean fragment size
     "-s",  "75", # std dev fragment size
     "-o",  output_prefix  # output prefix
     ]

logfile = output_prefix + ".art.log"
with open(logfile, 'w') as fh_log:
    print "Logging to " + logfile
    print "About to execute " + " ".join(ART_CMD)
    subprocess.check_call(ART_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)



# ART creates fastq files with suffix 1.fq.  We want .1.fq suffix.  Rename files.
shutil.move(output_prefix + "1.fq", output_prefix + ".1.fq")
shutil.move(output_prefix + "2.fq", output_prefix + ".2.fq")


# ART creates an error free sam.  Generate error free fq from error free sam
PICARD_CMD = ["java", "-jar", PICARD_BIN_DIR + os.sep + "SamToFastq.jar",
              "INPUT=" + output_prefix + "_errFree.sam",
              "FASTQ=" + output_prefix + ".errFree.1.fq",
              "SECOND_END_FASTQ=" + output_prefix + ".errFree.2.fq"]
picard_logfile = output_prefix + ".picard.log"
with open(picard_logfile, 'w') as fh_log:
    print "Logging to " + logfile
    print "About to execute " + " ".join(PICARD_CMD)
    subprocess.check_call(PICARD_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)


# Align reads against population consensus
bowtie_log = output_prefix + ".bowtie.log"
BOWTIE_BUILD_CMD = ["bowtie2-build",
             consensus_fasta,  # fasta to align against
             consensus_fasta.replace(".fasta", "")]  # prefix of db index

# Super high gap penalty, very sensitive local alignment
BOWTIE_CMD = ["bowtie2",
              "--very-sensitive-local",
              "-x", consensus_fasta.replace(".fasta", ""),  # reference index
              "-t",  # output timing
              "-1", output_prefix + ".1.fq",
              "-2", output_prefix + ".2.fq",
              "-S", output_prefix + ".consensus.bowtie.sam",  # samfile output
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
