"""
Helper module for renaming tree inner nodes.
"""
import Bio.Phylo as Phylo


def copy_node_names(src_treefile, dest_treefile, out_treefile):
    """
    For each clade (as determined by the tipnames) in src_treefile, finds the corresponding clade in rename_treefile,
    and renames the clade in treefile2.
    Writes out the modifed rename_treefile to out_treefile.

    src_treefile and rename_treefile must represent the same topology.  The nodes can be in different order,
    but the clades must contain the same tips in either tree.
    The only exception is that one tree can be rooted and the other tree unrooted.

    Does not rename the tips
    :return:
    """
    src_tree = Phylo.read(src_treefile, "newick")
    dest_tree = Phylo.read(dest_treefile, "newick")

    for src_clade in src_tree.find_clades(order="preorder"):
        if src_clade.is_terminal():
            continue

        dest_tips = []
        for src_tip in src_clade.get_terminals():
            dest_tip = dest_tree.find_clades(name=src_tip.name, terminal=True).next()
            dest_tips.append(dest_tip)

        dest_clade = dest_tree.common_ancestor(dest_tips)
        dest_clade.name = src_clade.name

    Phylo.write(dest_tree, out_treefile, "newick")