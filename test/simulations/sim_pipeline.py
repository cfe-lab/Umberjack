# The full pipeline for generating simulated population reads for unit testing.
# NB:   Needs picard v1.128 or later

import subprocess
import os
import logging
import sys
import ConfigParser
import indelible.indelible_handler as indelibler
import hyphy.hyphy_handler as hyphy_handler


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)


SECTION = "sim"

config_file = sys.argv[1]
config = ConfigParser.RawConfigParser()
config.read(config_file)

SEED = config.getint(SECTION, "SEED")
TREE_FILE_PREFIX = config.get(SECTION, "TREE_FILE_PREFIX")
NUM_CODON_SITES = config.getint(SECTION, "NUM_CODON_SITES")
NUM_INDIV = config.getint(SECTION, "NUM_INDIV")
OUTDIR = config.get(SECTION, "OUTDIR")
INDELIBLE_OUT_FILENAME_PREFIX = config.get(SECTION, "INDELIBLE_OUT_FILENAME_PREFIX")
INDELIBLE_BIN_DIR = config.get(SECTION, "INDELIBLE_BIN_DIR")
INTERVALS_CSV = config.get(SECTION, "INTERVALS_CSV")
SAMPLE_GENOMES_OUTDIR = config.get(SECTION, "SAMPLE_GENOMES_OUTDIR")
SAMPLE_GENOMES_OUT_FILENAME_PREFIX = config.get(SECTION, "SAMPLE_GENOMES_OUT_FILENAME_PREFIX")


# Generate Tree
asg_driver_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "asg_driver.py")
asg_driver_cmd = ["python", asg_driver_exe,
                  TREE_FILE_PREFIX,
                  str(NUM_INDIV),
                  str(SEED)]
LOGGER.debug("About to execute " + " ".join(asg_driver_cmd))
subprocess.check_call(asg_driver_cmd, env=os.environ)
LOGGER.debug("Finished execute ")
treefile = TREE_FILE_PREFIX + ".nwk"

# Relabel tree nodes to more manageable names.  Reformat tree so that indelible can handle it.
relabel_phylogeny_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "relabel_phylogeny.py")
relabel_phylogeny_cmd = ["python", relabel_phylogeny_exe,
                         treefile]
LOGGER.debug("About to execute " + " ".join(relabel_phylogeny_cmd))
subprocess.check_call(relabel_phylogeny_cmd, env=os.environ)
LOGGER.debug("Finished execute ")
renamed_treefile = TREE_FILE_PREFIX + ".rename.nwk"

# Use Indelible to create population sequences at different scaling factors (ie mutation rates)
os.environ["PATH"] += os.pathsep + INDELIBLE_BIN_DIR
batch_indelible_exe = os.path.abspath(os.path.dirname(__file__) + "/indelible/batch_indelible.py")
indelible_cmd = ["python", batch_indelible_exe,
                 renamed_treefile,  # full filepath to tree
                 INTERVALS_CSV,
                 str(SEED)]  # prefix of indelible population sequence output files
LOGGER.debug("About to execute " + " ".join(indelible_cmd))
subprocess.check_call(indelible_cmd, env=os.environ)
LOGGER.debug("Finished execute ")


# Create sample genome by concatenating slices of indelible alignments from different mutation rates.
sample_genomes_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "sample_genomes.py")
sample_genomes_cmd = ["python", sample_genomes_exe,
                      INTERVALS_CSV,  # full filepath of intervals CSV indicating indelible scaling factors and number of codons to take from each
                      SAMPLE_GENOMES_OUTDIR,  # full filepath of directory for sample_genomes.py output
                      SAMPLE_GENOMES_OUT_FILENAME_PREFIX]  # prefix of sample_genomes.py population sequence output files
LOGGER.debug("About to execute " + " ".join(sample_genomes_cmd))
subprocess.check_call(sample_genomes_cmd, env=os.environ)
LOGGER.debug("Finished execute ")


# Simulate MiSeq reads from the population genomes.
ART_BIN_DIR = config.get(SECTION, "ART_BIN_DIR")
ART_QUAL_PROFILE_TSV1 = config.get(SECTION, "ART_QUAL_PROFILE_TSV1")
ART_QUAL_PROFILE_TSV2 = config.get(SECTION, "ART_QUAL_PROFILE_TSV2")
ART_READ_OUTPUT_PREFIX =  config.get(SECTION, "ART_READ_OUTPUT_PREFIX")

PICARD_BIN_DIR = config.get(SECTION, "PICARD_BIN_DIR")
BOWTIE_OUT_DIR = config.get(SECTION, "BOWTIE_OUT_DIR")
PROCS = config.getint(SECTION, "PROCS")

sample_genomes_fasta = SAMPLE_GENOMES_OUTDIR + os.sep + SAMPLE_GENOMES_OUT_FILENAME_PREFIX + ".fasta"
sample_genomes_consensus_fasta = sample_genomes_fasta.replace(".fasta", ".consensus.fasta")

generate_reads_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "generate_reads.py")
generate_reads_cmd = ["python", generate_reads_exe,
                      ART_BIN_DIR,
                      ART_QUAL_PROFILE_TSV1,
                      ART_QUAL_PROFILE_TSV2,
                      sample_genomes_fasta,
                      sample_genomes_consensus_fasta,
                      ART_READ_OUTPUT_PREFIX,
                      PICARD_BIN_DIR,
                      BOWTIE_OUT_DIR,
                      str(PROCS),
                      str(SEED)]
LOGGER.debug("About to execute " + " ".join(generate_reads_cmd))
subprocess.check_call(generate_reads_cmd, env=os.environ)
LOGGER.debug("Finished execute ")




# Calculate HyPhy dN/dS against the full sample_genomes population fasta
# TODO:  do fasttree tree  since we lose the tree when we concatenate multiple populations at different scalings together


indelible_tree_txt = OUTDIR + os.sep + "50" + os.sep + "trees.txt"
indelible_out_tree_nwk = OUTDIR + os.sep + "50" + os.sep + INDELIBLE_OUT_FILENAME_PREFIX + ".50.nwk"
indelibler.write_tree(indelible_tree_txt, indelible_out_tree_nwk)
hyphy_handler.calc_dnds(codon_fasta_filename=sample_genomes_fasta,
                        tree_filename=indelible_out_tree_nwk,
                        pvalue=0.05,
                        threads=PROCS)


