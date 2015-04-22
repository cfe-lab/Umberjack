"""
Roots trees using java root to tip
"""

import os
import logging
import fasttree.fasttree_handler as fasttree
import Bio.Phylo as Phylo
import subprocess



LOGGER = logging.getLogger(__name__)


def make_rooted_tree(unrooted_treefile, threads=1):
    """
    Makes a tree and then roots the tree.  Also, names inner nodes and collapses zero-length branches into polytomies.
    :param unrooted_treefile:
    :return str:  path to rooted tree
    """

    # Expand polytomies so that we get a bifurcating tree.  otherwise ape.rtt() and richard's RootToTip will crash.
    bifurc_treefile = unrooted_treefile.replace(".nwk", ".2nary.nwk")
    LOGGER.debug("Removing polytomies in " + unrooted_treefile + " to " +  bifurc_treefile)
    if os.path.exists(bifurc_treefile) and os.path.getsize(bifurc_treefile):
        LOGGER.warn("Not regenerating " + bifurc_treefile)
    else:
        Rscript_wdir = os.path.dirname(os.path.realpath(__file__)) + "/../../R/timing"
        subprocess.check_call(["Rscript", Rscript_wdir + "/RemovePolytomies.r",
                               "--newick", "--outnewick",
                               unrooted_treefile, bifurc_treefile], env=os.environ)
        LOGGER.debug("Done removing polytomies in " + unrooted_treefile)


    # Custom RootToTip for rooting tree  (modified from the one from BEAST)
    rooted_treefile = unrooted_treefile.replace(".nwk", ".rooted.nxs")
    time_rooted_treefile = unrooted_treefile.replace(".nwk", ".rooted.time.nxs")
    node_dates_csv = unrooted_treefile.replace(".nwk", ".nodedates.csv")

    LOGGER.debug("Rooting tree to " + rooted_treefile + " and " + time_rooted_treefile)
    if os.path.exists(rooted_treefile) and os.path.getsize(rooted_treefile) and os.path.exists(time_rooted_treefile) and os.path.getsize(time_rooted_treefile):
        LOGGER.warn("Not regenerating rooted trees " + rooted_treefile + " and " + time_rooted_treefile)
    else:
        subprocess.check_call(["java", "-jar", os.path.dirname(os.path.realpath(__file__)) + os.sep + "RLRootToTip.jar",
                               "-writetree", rooted_treefile,
                               "-timetree", time_rooted_treefile,
                               "-nodedates", node_dates_csv,
                               "-newick",
                               bifurc_treefile],
                              env=os.environ)
        LOGGER.debug("Done rooting tree " + rooted_treefile)


    # If we leave the tree as binary tree with polytomies expanded by using zero-length branches,
    # when hyphy scales the tree branch lengths from nucleotide substitutions/site to codon substitutions/site,
    # it can set the zero-length branches to non-zero length codon branches, albeit with very very very small length.
    # In order to avoid this, we convert the binary tree back to polytomies before sending to hyphy to reconstruct ancestors.
    named_rooted_treefile = rooted_treefile.replace(".nxs", ".named.nwk")
    named_rooted_time_treefile = time_rooted_treefile.replace(".nxs", ".named.nwk")
    LOGGER.debug("Collapsing 0-length branches and adding inner node names to " + named_rooted_treefile + " and " + named_rooted_time_treefile)
    if (os.path.exists(named_rooted_treefile) and os.path.getsize(named_rooted_treefile) and
            os.path.exists(named_rooted_time_treefile) and os.path.getsize(named_rooted_time_treefile)):
        LOGGER.warn("Not regenerating " + named_rooted_treefile + " and " + named_rooted_time_treefile)
    else:
        rooted_tree = Phylo.read(rooted_treefile, "nexus")
        for cladenum, clade in enumerate(rooted_tree.find_clades(target=None, terminal=False, order="preorder")):
            if not clade.name:
                clade.name = "Clade" + str(cladenum)

        # TODO:  check if time scaled tree and nucleotide scaled tree have same topology
        rooted_time_tree = Phylo.read(time_rooted_treefile, "nexus")
        for cladenum, clade in enumerate(rooted_time_tree.find_clades(target=None, terminal=False, order="preorder")):
            if not clade.name:
                clade.name = "Clade" + str(cladenum)

        # Keep track of which clades we delete in nuc sub tree, since the timed tree can have zero-length branches where nuc tree doesn't
        # del_clades = list(rooted_tree.find_clades(lambda c: c.branch_length  == 0, terminal=False))
        # del_tips = []
        rooted_tree.collapse_all(lambda c: c.branch_length  == 0)  # Collapse all zero-length branches into polytomies
        for del_tip in rooted_tree.find_clades(lambda c: c.branch_length  == 0, terminal=True):
            if del_tip.branch_length == 0:
                parent_path = rooted_tree.get_path(del_tip)
                if len(parent_path) >=2:  # tip is not child of root
                    parent = parent_path[-2]
                    rooted_tree.collapse(parent)
                    # del_tips.extend([del_tip])


        Phylo.write(rooted_tree, named_rooted_treefile, "newick")
        Phylo.write(rooted_time_tree, named_rooted_time_treefile, "newick")
        LOGGER.debug("Done collapsing 0-length branches and adding inner node names to " + named_rooted_treefile + " and " + named_rooted_time_treefile)
        # Default biopython 1.58 with Ubuntu 12.04 doesn't write inner nodes names.  Need >biopython 1.59
        # If we output in biopython 1.65 nexus format, it uses str(clade) to get taxa labels, which cuts off the names and appends "...",
        #   which screws up hyphy parsing.
    return named_rooted_treefile
