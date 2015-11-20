"""
Invokes asg_driver.py to create an phylogeny following an ancestral selection graph.

Remove massive tip and internal node labels from asg_driver.py
simulated phylogeny.  INDELible doesn't like internal node labels.  It also doesn't like singeton root nodes.
Make tree suitable for INDELible input.

If we need to make variations of the same tree to simulate recombination,
then for B breakpoints, make B+1 trees.
At each breakpoint, the population's genomes will follow a random different tree.

"""
import sys
from Bio import Phylo
from cStringIO import StringIO
import random
import config.settings
import logging
from argparse import ArgumentParser
import os
import subprocess
import math
import recombination


LOGGER = logging.getLogger(__name__)

MAX_SEED = math.pow(2, 32)-1  # internally asg_driver.py used numpy random generator which can only take up to 32bit seeds


def make_tree(treefile_prefix, total_tips, seed, selection_rate, generations):
    """
    Invoke asg_driver to make tree
    :param str treefile_prefix:  filepath to treefile without the .nwk suffix
    :param int total_tips:  total tips in the tree
    :param int seed:  random seed
    :param float selection_rate:  rate that dominant allele is selected  (assuming 2 allele model in ancestral selection graph.  [0, 1]
    :param int generations:  total generations in which new individuals are introduced into the population.  Must be > 0
    :return str:  treefile name
    """
    treefile = treefile_prefix + ".nwk"

    if os.path.exists(treefile) and os.path.exists(treefile):
        LOGGER.warn("Not regenerating " + treefile)
    else:
        asg_driver_exe = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "asg_driver.py")
        asg_driver_cmd = ["python", asg_driver_exe,
                      "-f", treefile_prefix,
                      "-n", str(total_tips),
                      "-r", str(seed),
                      "-s", str(selection_rate),
                      "-g", str(generations)
                      ]

        LOGGER.debug("About to execute " + " ".join(asg_driver_cmd))
        subprocess.check_call(asg_driver_cmd, env=os.environ)
        LOGGER.debug("Finished execute " + " ".join(asg_driver_cmd))

    return treefile




def format_tree(format_treefile, treefile=None, tree=None):
    """
    Formats the treefiles so that they are compatible with INDELible.

    Renames the tip names to otu<tipnum> if they haven't been renamed already.  Removes inner node names.
    Collapses singleton root nodes.


    :param str format_treefile:  full filepath to output newick tree formatted for INDELible compatibility
    :param str treefile:  full filepath to input newick tree created by asg_driver.py.  Must be specified if tree unspecified.
    :param Bio.Phylo.Tree  tree:  tree instance if input tree filepath not given
    :return Bio.Phylo.Tree:  revised tree instance
    """

    if not treefile and not tree:
        raise ValueError("Must specify either treefile or tree")

    if tree:
        t = tree
    else:
        t = Phylo.read(treefile, 'newick')

    for i, node in enumerate(t.get_nonterminals()):
        node.name = "" #  indelible does not like named inner nodes

    for i, tip in enumerate(t.get_terminals()):
        if not tip.name.startswith("otu"):
            tip.name = 'otu'+str(i+1)

    # Indelible wants beginning of newick to start with (, end of newick to end with );
    # I.e.  (....);
    # It also doesn't like singleton root nodes.
    # R doesn't like singleton root nodes too.
    #
    # BioPython will write out root branch as singleton node no matter what.
    # EG.  ((...)Clade0:4, (...)Clade1:3)Root:2;
    # Collapse the root singleton node
    # EG>  ((...)Clade0:4, (...)Clade1:3);

    tree_strio = StringIO()
    Phylo.write(t, tree_strio, format="newick", format_branch_length='%1.9f')
    tree_strio.flush()
    tree_str = tree_strio.getvalue()
    tree_strio.close()
    # Find index of last branch length separator ":", this is the root branch length separator
    root_branch_index = tree_str.rindex(":")
    tree_str = tree_str[:root_branch_index] + ";"

    with open(format_treefile, "w") as fh_out:
        fh_out.write(tree_str)

    return t


def make_format_trees(treefile_prefixes, total_tips, selection_rate, generations, seed=None):
    """
    For each treefile, generate a separate tree.
    Format each tree so that it plays nice with INDELible.
    Each tree will be the same topology and total branch length as the previous tree except for a random prune-regraph of a random clade.

    :param list treefile_prefixes: comma separated list of treefile prefixes to write to.  Appends .rename.nwk to filepath of formatted trees.
        Appends .nwk to filepaths of unformatted trees.
    :param int total_tips:  number of tips in each tree
    :param int seed:  random seed
    :return list:  list of filepaths to each formatted breakpoint newick tree
    """

    # Use a custom random.Random instance when generating breakpoint trees.
    # so the the tree generation is not affected by calls to random for other purposes.
    # We seed the custom tree_randomizer with a known seed so that we can deterministicaly create the same
    # set of trees each time we run this script.
    tree_randomizer = random.Random(seed)

    # Ensure that each new tree is different from the last by seeding it with a new random int.
    # As long as the tree_randomizer instance is initialized with the same seed, we will always get the same set of new random ints
    # and thus the same set of new trees.
    treeseed = tree_randomizer.randint(0, MAX_SEED)
    unformat_treefile = make_tree(treefile_prefix=treefile_prefixes[0], total_tips=total_tips, seed=treeseed,
                                  selection_rate=selection_rate, generations=generations)

    # format newick files for INDELible compatibility
    format_treefile = unformat_treefile.replace(".nwk", ".rename.nwk")
    format_tree(treefile=unformat_treefile, format_treefile=format_treefile)


    # for remaining trees, randomly prune and regraph a clade from the previous tree
    prev_treefile = format_treefile
    for treefile_prefix in treefile_prefixes[1:]:
        treeseed = tree_randomizer.randint(0, MAX_SEED)

        graft_clade, recombo_tree = recombination.make_recombo_tree(in_treefile=prev_treefile, out_treefile=None, seed=treeseed)

        recombo_treefile = treefile_prefix + ".rename.nwk"
        format_tree(tree=recombo_tree, format_treefile=recombo_treefile)

        prev_treefile = recombo_treefile




if __name__ == "__main__":

    config.settings.setup_logging()  # Use default logging format from Umberjack settings

    # Parse commandline arguments
    parser = ArgumentParser()
    parser.add_argument("-f", help="Comma separated list of filename prefixes for the output newick tree files " +
                                   " (including filpath to directories) " +
                                   "Generates a new tree for each output file. "
                                   "For the unformatted trees, adds .nwk suffix.  For the formatted trees, adds .rename.nwk suffix")
    parser.add_argument("-t", help="Number of tips in each tree", type=int)
    parser.add_argument("-r", help="Integer to seed randomization of tree topology and branch length. " +
                                   "Each tree generated will be different from the previous, "
                                   "but the sequence of trees generated will be deterministic " +
                                   "between multiple processes running this script as long " +
                                   "as the seed used by each process is the same.", type=int, default=None)
    parser.add_argument("-s", help="selection rate used in ancestral selection graph to create tree.  Must be in [0, 1]", type=float)
    parser.add_argument("-g", help="generations used in ancestral selection graph to create tree.  Must be > 0", type=int)


    args = parser.parse_args()
    if not args:
        parser.print_usage()
        sys.exit()

    treefiles = args.f.split(",")

    make_format_trees(treefile_prefixes=treefiles, total_tips=args.t, seed=args.r,
                      selection_rate=args.s, generations=args.g)








