"""
Remove massive tip and internal node labels from Richard Liang's
simulated phylogeny.
While we're at it, we need to rescale the coalescence times 
(branch lengths) to confirm to a logistic growth model because 
the earliest coalescence events are way too deep.
"""
import sys
import math
from Bio import Phylo
from cStringIO import StringIO

in_tree_file = sys.argv[1]

t = Phylo.read(in_tree_file, 'newick')


nodes = t.get_nonterminals()
for i, node in enumerate(nodes):
    node.name = "" #  "N" + str(i)  # indelible does not like named inner nodes

tips = t.get_terminals()
for i, tip in enumerate(tips):
    tip.name = 'otu'+str(i+1)



#t.collapse(t.clade)

# # dictionary of Clade to depth (tree height)
# # depths = t.depths()
#
# #max_height = max(depths.values())
#
# # apply a shrink factor that increases with depth
# # in a logistic manner, mimicking logistic growth
# branches = nodes + tips
#
# kcap = 1000.
# n0 = 30.
# r = 20.
#
# global_scaling_factor = 10000.
#
# # # adjust coalescent times with logistic growth factor
# for i, branch in enumerate(branches):
#     cur = branch.branch_length
#     #depth = 10.* (max_height - depths[branch]) / max_height # scaled to range from 0 to 10.
#     #logistic = kcap * n0 * math.exp(r * depth) / (kcap + n0 * (math.exp(r * depth) - 1))
#     #adj = cur / logistic
#     #branch.branch_length = adj / global_scaling_factor
#     branch.branch_length = cur / global_scaling_factor
#     #if i % 500 == 0: print i, cur, depth, logistic, adj


# branches = nodes + tips
# t.clade.branch_length = 0.0
# total_branch_len = t.total_branch_length()
# global_scaling_factor = 100.0
# for i, branch in enumerate(branches):
#     branch_scaling_factor = global_scaling_factor/total_branch_len
#     branch.branch_length = branch.branch_length * branch_scaling_factor


# Indelible wants beginning of newick to start with (, end of newick to end with );
# It also doesn't like singleton root nodes.
# R doesn't like singleton root nodes too.

# BioPython will write out root branch as singleton node no matter what.
# EG.  ((...)Clade0:4, (...)Clade1:3)Root:2;
# Collapse the root single node
# EG>  ((...)Clade0:4, (...)Clade1:3);
tree_strio = StringIO()
Phylo.write(t, tree_strio, format="newick", format_branch_length='%1.9f')
tree_strio.flush()
tree_str = tree_strio.getvalue()
tree_strio.close()
# Find index of last branch length separator ":", this is the root branch length separator
root_branch_index = tree_str.rindex(":")
tree_str = tree_str[:root_branch_index] + ";"



relabled_nwk = in_tree_file.replace(".nwk", ".rename.nwk")
with open(relabled_nwk, "w") as fh_out:
    fh_out.write(tree_str)

#Phylo.write(t, file=in_tree_file.replace(".nwk", ".rename.nwk"), format='newick', format_branch_length='%1.9f')
