import logging
import subprocess
import os
import config.settings as settings


LOGGER = logging.getLogger(__name__)


ENV_OMP_NUM_THREADS = 'OMP_NUM_THREADS'
GTRRATES_LINE_START = "GTRRates"



def make_tree(fasta_fname, threads=1, fastree_exe=settings.DEFAULT_FASTTREEMP_EXE, debug=False, custom_flags=None):
    """
    Creates a phylogenetic tree from the fasta.
    Tree file written to same directory as fasta and has same name as fasta but with ".nwk" suffix.
    FastTree logfile written to same directory as fasta and has same name as fasta but with ".fasttree.log" suffix.
    :param str fasta_fname:  path to multiple sequence aligned fasta file.  File should have .fasta suffix.
    :param int threads: number of FastTreeMP threads
    :param str fastree_exe:  path to FastTreeMP or FastTree executable.  If empty, then uses FastTreeMP from PATH environment variable.
    :param bool debug:  whether to output verbose fasttree stdout/stderr to file.  File will have same name as fasta but with ".fasttree.stdouterr.txt" suffix.
    :param list custom_flags:  list of strings for each custom FastTreeMP commandline argument
    :return : file path to tree
    :rtype :  str
    """
    fasta_fname_prefix = os.path.splitext(fasta_fname)[0]  # Remove the .fasta suffix
    fastree_logfilename = fasta_fname_prefix + ".fasttree.log"
    fastree_treefilename = fasta_fname_prefix + ".nwk"
    fasttree_stdouterr_filename = fasta_fname_prefix + ".fasttree.stdouterr.txt"
    LOGGER.debug("Start Fasttree  " + fastree_treefilename)

    if not os.path.exists(fastree_treefilename) or os.path.getsize(fastree_treefilename) <= 0:
        os.environ[ENV_OMP_NUM_THREADS] = str(threads)
        fasttree_cmd = [fastree_exe,
                        '-gtr',  # general time reversible nucleotide substitution model
                        '-nt',   # nucleotides
                        '-gamma',  # sites vary with 20 category gamma distro
                        '-nosupport',  # do not output support values in tree
                        '-log', fastree_logfilename,  # fast tree log
                        '-out', fastree_treefilename,  # tree
        ]

        if custom_flags:
            fasttree_cmd = fasttree_cmd + custom_flags

        if debug:
            with open(fasttree_stdouterr_filename, 'w') as fasttree_stdouterr_fh:
                fasttree_cmd = fasttree_cmd + [fasta_fname]  # multiple sequence aligned fasta
                subprocess.check_call(fasttree_cmd,
                                  stdout=fasttree_stdouterr_fh, stderr=fasttree_stdouterr_fh, shell=False,
                                  env=os.environ)
        else:
            with open(os.devnull, 'wb') as fh_devnull:
                fasttree_cmd = fasttree_cmd + [
                    '-quiet',  # suppress reporting information
                    '-nopr',  # suppress progress indicator
                    fasta_fname]  # multiple sequence aligned fasta
                subprocess.check_call(fasttree_cmd,  # multiple sequence aligned fasta
                                      stdout=fh_devnull, stderr=fh_devnull, shell=False,
                                      env=os.environ)
        LOGGER.debug("Done Fasttree " + fastree_treefilename)
    else:
        LOGGER.debug("Found existing Fasttree " + fastree_treefilename + ". Not regenerating")

    return fastree_treefilename


def make_tree_repro(fasta_fname, intree_fname, fastree_exe=settings.DEFAULT_FASTREE_EXE):
    """
    Creates a phylogenetic tree from the fasta.  Resulting tree will maintain the same topology as input tree,
    but branch lengths will be optimized.

    NB:  FastTree2.1 sometimes changes the topology if there are duplicate sequences.

    Tree file written to same directory as fasta and has same name as fasta but with ".nwk" suffix.
    FastTree logfile written to same directory as fasta and has same name as fasta but with ".fasttree.log" suffix.
    Verbose FastTree stdout/stderr written to same directory as fasta but with ".fasttree.stdouterr.txt" suffix.

    :param str fasta_fname:  multiple sequence aligned fasta file
    :param str intree_fname:  full filepath to input tree for which the topology should not change
    :param str fastree_exe: full filepath to FastTree executable.  If empty, uses FastTree from PATH env var
    """
    fasta_fname_prefix = os.path.splitext(fasta_fname)[0]  # Remove the .fasta suffix
    fastree_logfilename = fasta_fname_prefix + ".fasttree.log"
    fastree_treefilename = fasta_fname_prefix + ".nwk"
    fasttree_stdouterr_filename = fasta_fname_prefix + ".fasttree.stdouterr.txt"
    LOGGER.debug("Start Fasttree  " + fastree_treefilename)

    if not os.path.exists(fastree_treefilename) or os.path.getsize(fastree_treefilename) <= 0:
        with open(fasttree_stdouterr_filename, 'w') as fasttree_stdouterr_fh:
            subprocess.check_call([fastree_exe,
                                   '-gtr',  # general time reversible nucleotide substitution model
                                   '-nt',   # nucleotides
                                   '-gamma',  # sites vary with 20 category gamma distro
                                   '-nosupport',  # do not output support values in tree
                                   "-nome", # no minimum evolution nearest neighbour exchange or subtree prune regraph
                                   "-mllen",  # find max likelihood branch lengths
                                   "-intree", intree_fname,  #  input tree newick file.  Topology of resulting tree will not change.
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
    :param str fastree_logfilename:  path to fasttree logfile
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
                str_rates = line.rstrip().split()[1:]
                (AC, AG, AT, CG, CT, GT) = (float(x) for x in str_rates)
                return AC, AG, AT, CG, CT, GT





