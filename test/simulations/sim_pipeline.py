"""
The full pipeline for generating simulated population reads for unit testing.
Usage:  python sim_pipeline.py [config file]
"""

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


def get_path_str(path, pardir):
    """
    If absolute path, then returns the path as is.
    If relative path, then returns absolute path of concatenated pardir/path
    :param str path:  absolute or relative file or directory path
    :param str pardir: parent directory to concatenate to path if path is relative directory
    :return str: absolute resolved path
    """
    if not os.path.isabs(path):
        return os.path.join(pardir, path)
    else:
        return path



SECTION = "sim"

config_file = sys.argv[1]
config = ConfigParser.RawConfigParser()
config.read(config_file)

OUTDIR = os.path.dirname(config_file)  # Output directory for simulated data



# Generate Tree
SEED = config.getint(SECTION, "SEED")
FILENAME_PREFIX = config.get(SECTION, "FILENAME_PREFIX")
NUM_CODON_SITES = config.getint(SECTION, "NUM_CODON_SITES")
NUM_INDIV = config.getint(SECTION, "NUM_INDIV")

treefile = OUTDIR + os.sep + FILENAME_PREFIX + ".nwk"
renamed_treefile = OUTDIR + os.sep + FILENAME_PREFIX + ".rename.nwk"
if os.path.exists(treefile) and os.path.getsize(treefile) and os.path.exists(renamed_treefile) and os.path.getsize(renamed_treefile):
    LOGGER.warn("Not regenerating trees {} and {}".format(treefile, renamed_treefile) )
else:
    asg_driver_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "asg_driver.py")
    asg_driver_cmd = ["python", asg_driver_exe,
                      OUTDIR + os.sep + FILENAME_PREFIX,
                      str(NUM_INDIV),
                      str(SEED)]
    LOGGER.debug("About to execute " + " ".join(asg_driver_cmd))
    subprocess.check_call(asg_driver_cmd, env=os.environ)
    LOGGER.debug("Finished execute ")


    # Relabel tree nodes to more manageable names.  Reformat tree so that indelible can handle it.
    relabel_phylogeny_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "relabel_phylogeny.py")
    relabel_phylogeny_cmd = ["python", relabel_phylogeny_exe,
                             treefile]
    LOGGER.debug("About to execute " + " ".join(relabel_phylogeny_cmd))
    subprocess.check_call(relabel_phylogeny_cmd, env=os.environ)
    LOGGER.debug("Finished execute ")



# Use Indelible to create population sequences at different scaling factors (ie mutation rates)
INDELIBLE_BIN_DIR = get_path_str(config.get(SECTION, "INDELIBLE_BIN_DIR"), OUTDIR)
INDELIBLE_SCALING_RATES = config.get(SECTION, "INDELIBLE_SCALING_RATES")

batch_indelible_exe = os.path.abspath(os.path.dirname(__file__) + "/indelible/batch_indelible.py")
indelible_cmd = ["python", batch_indelible_exe,
                 renamed_treefile,  # full filepath to tree
                 INDELIBLE_SCALING_RATES,
                 str(SEED),  # random seed
                 str(NUM_CODON_SITES), # number of codon sites in genome
                 OUTDIR,  # indelible output file directory
                 FILENAME_PREFIX,  # Indelible output filename prefix
                 INDELIBLE_BIN_DIR]  # indelible bin dir
LOGGER.debug("About to execute " + " ".join(indelible_cmd))
subprocess.check_call(indelible_cmd, env=os.environ)
LOGGER.debug("Finished execute ")


# Create sample genome by concatenating slices of indelible alignments from different mutation rates.
sample_genomes_fasta = OUTDIR + os.sep + "mixed" + os.sep + FILENAME_PREFIX + ".mixed.fasta"
sample_genomes_consensus_fasta = sample_genomes_fasta.replace(".fasta", ".consensus.fasta")
if (os.path.exists(sample_genomes_fasta) and os.path.getsize(sample_genomes_fasta) and
        os.path.exists(sample_genomes_consensus_fasta) and os.path.getsize(sample_genomes_consensus_fasta)):
    LOGGER.warn("Not regenerating combined sample genome fastas {} and {} ".format(sample_genomes_fasta, sample_genomes_consensus_fasta))
else:
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
generate_reads_exe = os.path.abspath(os.path.dirname(__file__) + os.sep + "generate_reads.py")
generate_reads_cmd = ["python", generate_reads_exe,
                      config_file]
LOGGER.debug("About to execute " + " ".join(generate_reads_cmd))
subprocess.check_call(generate_reads_cmd, env=os.environ)
LOGGER.debug("Finished execute ")

# For the sample_genomes populations, we lose the true tree branch lengths when we concatenate multiple populations at different scalings together.
# Get FastTree to approximate tree for concatenated population sequences.
FASTTREE_EXE = get_path_str(config.get(SECTION, "FASTTREE_EXE"), OUTDIR)
sample_genomes_tree_fname = fasttree_handler.make_tree_repro(fasta_fname=sample_genomes_fasta, intree_fname=renamed_treefile,
                                                             fastree_exe=FASTTREE_EXE)


# Calculate HyPhy dN/dS for the full sample_genomes population fasta
PROCS = config.getint(SECTION, "PROCS")
HYPHY_EXE = get_path_str(config.get(SECTION, "HYPHY_EXE"), OUTDIR)
HYPHY_BASEPATH = get_path_str(config.get(SECTION, "HYPHY_BASEPATH"), OUTDIR)
hyphy_handler.calc_dnds(codon_fasta_filename=sample_genomes_fasta, tree_filename=sample_genomes_tree_fname,
                        hyphy_exe=HYPHY_EXE, hyphy_basedir=HYPHY_BASEPATH, threads=PROCS)


