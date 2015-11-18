"""
Roots trees using java root to tip
"""

import os
import logging
import fasttree.fasttree_handler as fasttree
import Bio.Phylo as Phylo
import subprocess
import Utility


LOGGER = logging.getLogger(__name__)


def make_rooted_tree(unrooted_treefile, threads=1):
    """
    Makes a tree and then roots the tree.  Also, names inner nodes and collapses zero-length branches into polytomies.
    :param unrooted_treefile:
    :return str:  path to rooted tree
    """
    # Expand polytomies so that we get a bifurcating tree.  otherwise ape.rtt() and RootToTip will crash.
    file_prefix = os.path.splitext(unrooted_treefile)[0]  # Remove the file suffix
    bifurc_treefile = file_prefix + ".binary.nwk"
    LOGGER.debug("Removing polytomies in " + unrooted_treefile + " to " +  bifurc_treefile)
    if os.path.exists(bifurc_treefile) and os.path.getsize(bifurc_treefile):
        LOGGER.warn("Not regenerating " + bifurc_treefile)
    else:
        Rscript_wdir = os.path.dirname(os.path.realpath(__file__)) + "/../R/timing"
        subprocess.check_call(["Rscript", Rscript_wdir + "/RemovePolytomies.r",
                               "--newick", "--outnewick",
                               unrooted_treefile, bifurc_treefile], env=os.environ)
        LOGGER.debug("Done removing polytomies in " + unrooted_treefile)


    # Custom RootToTip for rooting tree  (modified from the one from BEAST)
    rooted_treefile = file_prefix + ".rooted.nxs"
    time_rooted_treefile = file_prefix + ".rooted.time.nxs"
    node_dates_csv = file_prefix + ".nodedates.csv"
    root_results_csv = file_prefix + ".root.result.csv"
    LOGGER.debug("Rooting tree to " + rooted_treefile + " and " + time_rooted_treefile)
    if os.path.exists(rooted_treefile) and os.path.getsize(rooted_treefile) and os.path.exists(time_rooted_treefile) and os.path.getsize(time_rooted_treefile):
        LOGGER.warn("Not regenerating rooted trees " + rooted_treefile + " and " + time_rooted_treefile)
    else:
        subprocess.check_call(["java", "-jar", os.path.dirname(os.path.realpath(__file__)) + os.sep + "bin" + os.sep + "RLRootToTip.jar",
                               "-writetree", rooted_treefile,
                               "-timetree", time_rooted_treefile,
                               "-nodedates", node_dates_csv,
                               "-newick",
                               bifurc_treefile,
                               root_results_csv],
                              env=os.environ)
        LOGGER.debug("Done rooting tree " + rooted_treefile)


    # If we leave the tree as binary tree with polytomies expanded by using zero-length branches,
    # when hyphy scales the tree branch lengths from nucleotide substitutions/site to codon substitutions/site,
    # it can set the zero-length branches to non-zero length codon branches, albeit with very very very small length.
    # In order to avoid this, we convert the binary tree back to polytomies before sending to hyphy to reconstruct ancestors.
    named_rooted_treefile = file_prefix + ".named.nwk"
    named_rooted_time_treefile = file_prefix + ".named.nwk"
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
        del_clades = list(rooted_tree.find_clades(lambda c: c.branch_length  == 0, terminal=False))
        del_tips = []
        rooted_tree.collapse_all(lambda c: c.branch_length  == 0)  # Collapse all zero-length inner node branches into polytomies


        # collapse all zero-length terminal branches into polytomies
        for del_tip in list(rooted_tree.find_clades(lambda c: c.branch_length  == 0, terminal=True)):
            if del_tip.branch_length == 0:
                parent_path = rooted_tree.get_path(del_tip)
                if len(parent_path) >=2:  # tip is not child of root
                    parent = parent_path[-2]
                    rooted_tree.collapse(parent)
                    del_tips.extend([del_tip])

        # Don't call rooted_time_tree.collapse_all() because the timed tree can have different number of zero-length branches than nuc sub tree
        for del_clade in del_clades:
            time_del_clades = list(rooted_time_tree.find_clades(name=del_clade.name, terminal=False))
            if len(time_del_clades) > 1:
                raise ValueError("Timed tree has multiple nodes with name " + del_clade.name)

            if time_del_clades[0] != rooted_time_tree.root:  # collapse_all() never collapses root.  Since rooted_tree has its root intact, so should rooted_time_tree.
                rooted_time_tree.collapse(time_del_clades[0])

        for del_tip in del_tips:
            time_del_tips = list(rooted_time_tree.find_clades(name=del_tip.name, terminal=True))
            if len(time_del_tips) > 1:
                raise ValueError("Timed tree has multiple tips with name " + del_tip.name)
            parent_path = rooted_time_tree.get_path(time_del_tips[0])
            if len(parent_path) >=2:  # tip is not child of root
                parent = parent_path[-2]
                rooted_time_tree.collapse(parent)

        Phylo.write(rooted_tree, named_rooted_treefile, "newick")
        Phylo.write(rooted_time_tree, named_rooted_time_treefile, "newick")
        LOGGER.debug("Done collapsing 0-length branches and adding inner node names to " + named_rooted_treefile + " and " + named_rooted_time_treefile)
        # Default biopython 1.58 with Ubuntu 12.04 doesn't write inner nodes names.  Need >biopython 1.59
        # If we output in biopython 1.65 nexus format, it uses str(clade) to get taxa labels, which cuts off the names and appends "...",
        #   which screws up hyphy parsing.
    return named_rooted_treefile


def is_timeable(msa_fasta):
    """
    Checks if the are multiple timepoints.
    Reads should have "_<time from baseline>" in their names to specify the timepoint from which they came.
    :param str msa_fasta:  filepath to multiple-sequence aligned fasta
    :return bool:  whether the fasta contains reads from multiple timepoints
    :raises ValueError:  if the file is missing or does not follow the required timepoint readname format.
    """
    if not os.path.exists(msa_fasta):
        raise ValueError("Missing file " + msa_fasta)

    elif not os.path.getsize(msa_fasta):
        return False

    # Assume that reads from the same timepoint are aggregated together in file.
    # Read in first and last header so that we don't have to iterate through entire file to find different timepoints.
    # If that fails, then read header by header.
    first_header = Utility.get_first_header(msa_fasta)
    last_header = Utility.get_last_header(msa_fasta)
    if not first_header or not last_header:
        raise ValueError("Invalid multiple-sequence alignment fasta format " + msa_fasta)

    try:
        first_timepoint = float(first_header.split("_")[-1])
        last_timepoint = float(last_header.split("_")[-1])
    except ValueError:
        raise ValueError("Either first or last read do not follow required name format <name>_<time since baseline> in " + msa_fasta)

    if first_timepoint != last_timepoint:
        return True

    prev_timepoint = None
    for i, header in enumerate(Utility.get_fasta_headers(msa_fasta)):
        try:
            timepoint = float(header.split("_")[-1])
        except ValueError:
            raise ValueError("Read  " + str(i+1)  + " does not follow required name format <name>_<time since baseline> in " + msa_fasta)

        if prev_timepoint is not None and timepoint != prev_timepoint:
            return True
        prev_timepoint = timepoint

    return False