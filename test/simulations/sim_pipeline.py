# The full pipeline for generating simulated population reads for unit testing.

#
# * Create population phylogeny:
#
# ```bash
# python asg_driver.py  [newick tree filename prefix]
# ```
#
# * Create multiple sets of population DNA sequence, with each set using different branch scaling factor.
# 	* Install [indelible](http://abacus.gene.ucl.ac.uk/software/indelible/) and put the install folder in your PATH environment variable
# 	* You can generate multiple populations.  The amino acid substitution rates will follow a fixed 61-category gamma distribution.  However, you can set the branch scaling factor for each population.
#
# ```bash
# python batch_indelible.py [newick tree full filename] [comma-delimited list of scaling factors]
# ```
# 	* For each population, genome DNA sequences will be generated in scaling_<scaling factor>_TRUE.fas.
# 	* For each population, the amino acid substitution rate category at each site will be generated in scaling_<scaling factor>_RATES.txt
#
#
# * Create final population DNA sequence by concatenating the sections of previous DNA sequence from different scaling factors.
#
# ```bash
# python sample_genomes.py
# ```
#
# * Simulate miseq reads for final population with (ART)[http://www.niehs.nih.gov/research/resources/software/biostatistics/art/], a read simulator tool

import subprocess
import os
import logging
import sys
import ConfigParser
import shutil

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

TREE_FILE = config.get(SECTION, "TREE_FILE")
SCALING_FACTOR = config.getfloat(SECTION, "SCALING_FACTOR")
NUM_CODON_SITES = config.getint(SECTION, "NUM_CODON_SITES")
OUTDIR = config.get(SECTION, "OUTDIR")
OUTPUT_FILENAME_PREFIX = config.get(SECTION, "OUTPUT_FILENAME_PREFIX")
INDELIBLE_BIN_DIR = config.get(SECTION, "INDELIBLE_BIN_DIR")

renamed_treefile = TREE_FILE.replace(".nwk", ".rename.nwk")

# TODO:  actually rename the nodes properly
if not os.path.exists(renamed_treefile):
    shutil.copy(TREE_FILE, renamed_treefile)


os.environ["PATH"] += os.pathsep + INDELIBLE_BIN_DIR
scaled_popn_dir = OUTDIR + os.sep + str(SCALING_FACTOR)
batch_indelible_exe = os.path.abspath(os.path.dirname(__file__) + "/indelible/batch_indelible.py")
indelible_cmd = ["python", batch_indelible_exe,
                 renamed_treefile,  # full filepath to tree
                 str(SCALING_FACTOR),  # global scaling factor to stretch tree
                 str(NUM_CODON_SITES),  # number of codon sites in genome
                 scaled_popn_dir,  # output directory for indelible population sequences with mutation rate scaled by scaling factor
                 OUTPUT_FILENAME_PREFIX]  # prefix of indelible population sequence output files
LOGGER.debug("About to execute " + " ".join(indelible_cmd))
subprocess.check_call(indelible_cmd, env=os.environ)
LOGGER.debug("Finished execute ")
