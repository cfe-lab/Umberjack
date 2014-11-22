import logging
import subprocess
import os
import sys

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)

FASTTREE_EXE = "FastTreeMP"
ENV_OMP_NUM_THREADS = 'OMP_NUM_THREADS'
GTRRATES_LINE_START = "GTRRates"


def make_tree(fasta_fname, threads=1, fastree_exe=FASTTREE_EXE):
    """
    Creates a phylogenetic tree from the fasta
    :param str fasta_fname:  multiple sequence aligned fasta file
    :param int threads: number of threads allotted to processing this window  (only FastTree and HyPhy will be multithreaded)
    :param str fastree_exe: full filepath to FastTreeMP executable.  Default uses FastTreeMP from PATH env var
    """
    fasta_fname_prefix = os.path.splitext(fasta_fname)[0]  # Remove the .fasta suffix
    fastree_logfilename = fasta_fname_prefix + ".fasttree.log"
    fastree_treefilename = fasta_fname_prefix + ".tree"
    fasttree_stdouterr_filename = fasta_fname_prefix + ".fasttree.stdouterr.txt"
    LOGGER.debug("Start Fasttree  " + fastree_treefilename)

    if not os.path.exists(fastree_treefilename) or os.path.getsize(fastree_treefilename) <= 0:
        os.environ[ENV_OMP_NUM_THREADS] = str(threads)
        with open(fasttree_stdouterr_filename, 'w') as fasttree_stdouterr_fh:
            subprocess.check_call([fastree_exe,
                                   '-gtr',  # general time reversible nucleotide substitution model
                                   '-nt',   # nucleotides
                                   '-gamma',  # sites vary with 20 category gamma distro
                                   '-nosupport',  # do not output support values in tree
                                   '-log', fastree_logfilename,  # fast tree log
                                   '-out', fastree_treefilename,  # fast tree STDOUT, STDERR
                                   fasta_fname],  # multiple sequence aligned fasta
                                  stdout=fasttree_stdouterr_fh, stderr=fasttree_stdouterr_fh, shell=False,
                                  env=os.environ)
        LOGGER.debug("Done Fasttree " + fastree_treefilename)
    else:
        LOGGER.debug("Found existing Fasttree " + fastree_treefilename + ". Not regenerating")

    return fastree_treefilename


def extract_gtr_rates(fastree_logfilename):
    """
    Extracts the general time reversible nucleotide substution model rates from the fasttree log.
    Note:  FastTree gives rates with respect to GT
    :return tuple:  (AC, AG, AT, CG, CT, GT)
    """
    if not os.path.exists(fastree_logfilename) or not os.path.getsize(fastree_logfilename):
        LOGGER.error("FastTree Log " + fastree_logfilename + " is empty or does not exist")
        return None

    # parse the fasttree logfile for the GTRRates
    # EG)
    # GTRFreq	0.3654	0.1708	0.2532	0.2105
    # GTRRates	2.9049	15.7173	1.5149	0.2436	23.4565	1.0000
    with open(fastree_logfilename, 'rU') as fh:
        for line in fh:
            if line.startswith(GTRRATES_LINE_START):
                str_rates = line[len(GTRRATES_LINE_START):].lstrip().rstrip().split()
                (AC, AG, AT, CG, CT, GT) = (float(x) for x in str_rates)
                return AC, AG, AT, CG, CT, GT





