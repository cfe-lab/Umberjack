# The full pipeline for generating simulated population reads for unit testing.
# NB:   Needs picard v1.128 or later

import subprocess
import os
import logging
import sys
import ConfigParser
import hyphy.hyphy_handler as hyphy_handler
import fasttree.fasttree_handler as fasttree_handler
import config.settings as settings

settings.setup_logging()
LOGGER = logging.getLogger(__name__)



SECTION = "sim"

config_file = sys.argv[1]
config = ConfigParser.RawConfigParser()
config.read(config_file)

OUTDIR = os.path.dirname(config_file)  # Output directory for simulated data


SEED = config.getint(SECTION, "SEED")
FILENAME_PREFIX = config.get(SECTION, "FILENAME_PREFIX")
NUM_CODON_SITES = config.getint(SECTION, "NUM_CODON_SITES")
NUM_INDIV = config.getint(SECTION, "NUM_INDIV")


INDELIBLE_BIN_DIR = config.get(SECTION, "INDELIBLE_BIN_DIR")
INDELIBLE_SCALING_RATES = config.get(SECTION, "INDELIBLE_SCALING_RATES")

# Generate Tree
asg_driver_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "asg_driver.py")
asg_driver_cmd = ["python", asg_driver_exe,
                  OUTDIR + os.sep + FILENAME_PREFIX,
                  str(NUM_INDIV),
                  str(SEED)]
LOGGER.debug("About to execute " + " ".join(asg_driver_cmd))
subprocess.check_call(asg_driver_cmd, env=os.environ)
LOGGER.debug("Finished execute ")
treefile = OUTDIR + os.sep + FILENAME_PREFIX + ".nwk"

# Relabel tree nodes to more manageable names.  Reformat tree so that indelible can handle it.
relabel_phylogeny_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "relabel_phylogeny.py")
relabel_phylogeny_cmd = ["python", relabel_phylogeny_exe,
                         treefile]
LOGGER.debug("About to execute " + " ".join(relabel_phylogeny_cmd))
subprocess.check_call(relabel_phylogeny_cmd, env=os.environ)
LOGGER.debug("Finished execute ")
renamed_treefile = OUTDIR + os.sep + FILENAME_PREFIX + ".rename.nwk"


# Use Indelible to create population sequences at different scaling factors (ie mutation rates)
os.environ["PATH"] += os.pathsep + INDELIBLE_BIN_DIR
batch_indelible_exe = os.path.abspath(os.path.dirname(__file__) + "/indelible/batch_indelible.py")
indelible_cmd = ["python", batch_indelible_exe,
                 renamed_treefile,  # full filepath to tree
                 INDELIBLE_SCALING_RATES,
                 str(SEED),  # random seed
                 str(NUM_CODON_SITES), # number of codon sites in genome
                 OUTDIR,  # indelible output file directory
                 FILENAME_PREFIX]  # Indelible output filename prefix
LOGGER.debug("About to execute " + " ".join(indelible_cmd))
subprocess.check_call(indelible_cmd, env=os.environ)
LOGGER.debug("Finished execute ")


# Create sample genome by concatenating slices of indelible alignments from different mutation rates.
sample_genomes_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "sample_genomes.py")
sample_genomes_cmd = ["python", sample_genomes_exe,
                      INDELIBLE_SCALING_RATES,  #  comma delimited list of mutation scaling rates
                      OUTDIR + os.sep + "mixed",  # full filepath of directory for sample_genomes.py output
                      FILENAME_PREFIX + ".mixed", # prefix of sample_genomes.py population sequence output files
                      str(SEED), # random seed
                      str(NUM_CODON_SITES),  # number codon sites
                      OUTDIR,  # Indelible output directory
                      FILENAME_PREFIX]  # INDELible output filename prefix
LOGGER.debug("About to execute " + " ".join(sample_genomes_cmd))
subprocess.check_call(sample_genomes_cmd, env=os.environ)
LOGGER.debug("Finished execute ")


# Simulate MiSeq reads from the population genomes.
ART_BIN_DIR = config.get(SECTION, "ART_BIN_DIR")
ART_QUAL_PROFILE_TSV1 = config.get(SECTION, "ART_QUAL_PROFILE_TSV1")
ART_QUAL_PROFILE_TSV2 = config.get(SECTION, "ART_QUAL_PROFILE_TSV2")
ART_FOLD_COVER = config.getint(SECTION, "ART_FOLD_COVER")
ART_MEAN_INSERT = config.getint(SECTION, "ART_MEAN_INSERT")
ART_STDEV_INSERT = config.getint(SECTION, "ART_STDEV_INSERT")

PICARD_BIN_DIR = config.get(SECTION, "PICARD_BIN_DIR")
BWA_BIN_DIR = config.get(SECTION, "BWA_BIN_DIR")

PROCS = config.getint(SECTION, "PROCS")

sample_genomes_fasta = OUTDIR + os.sep + "mixed" + os.sep + FILENAME_PREFIX + ".mixed.fasta"
sample_genomes_consensus_fasta = sample_genomes_fasta.replace(".fasta", ".consensus.fasta")

art_reads_dir = OUTDIR + os.sep + "mixed" + os.sep + "reads"
art_reads_filename_prefix = FILENAME_PREFIX + ".mixed.reads"
generate_reads_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "generate_reads.py")
generate_reads_cmd = ["python", generate_reads_exe,
                      ART_BIN_DIR,
                      ART_QUAL_PROFILE_TSV1,
                      ART_QUAL_PROFILE_TSV2,
                      sample_genomes_fasta,
                      sample_genomes_consensus_fasta,
                      art_reads_dir + os.sep + art_reads_filename_prefix,  # dir and filename prefix of ART output
                      str(ART_FOLD_COVER),
                      str(ART_MEAN_INSERT),
                      str(ART_STDEV_INSERT),
                      PICARD_BIN_DIR,
                      BWA_BIN_DIR,
                      OUTDIR + os.sep + "mixed" + os.sep + "aln",  # BWA output dir
                      str(PROCS),
                      str(SEED),
                      OUTDIR + os.sep + "mixed" + os.sep + FILENAME_PREFIX + ".rates.csv"]  # Indelible mixed mutation rates csv
LOGGER.debug("About to execute " + " ".join(generate_reads_cmd))
subprocess.check_call(generate_reads_cmd, env=os.environ)
LOGGER.debug("Finished execute ")

# For the sample_genomes populations, we lose the true tree branch lengths when we concatenate multiple populations at different scalings together.
# Get FastTree to approximate tree for concatenated population sequences.
FASTTREE_EXE = config.get(SECTION, "FASTTREE_EXE")
sample_genomes_tree_fname = fasttree_handler.make_tree_repro(fasta_fname=sample_genomes_fasta, intree_fname=renamed_treefile,
                                                             fastree_exe=FASTTREE_EXE)


# Calculate HyPhy dN/dS for the full sample_genomes population fasta
HYPHY_EXE = config.get(SECTION, "HYPHY_EXE")
HYPHY_BASEPATH = config.get(SECTION, "HYPHY_BASEPATH")
hyphy_handler.calc_dnds(codon_fasta_filename=sample_genomes_fasta, tree_filename=sample_genomes_tree_fname,
                        hyphy_exe=HYPHY_EXE, hyphy_basedir=HYPHY_BASEPATH, threads=PROCS)


