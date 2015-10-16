# Checks tree stability from sequence subsampling.
# Do we get the same tree as the full population tree minus unsampled reads when we slide windows across the genome?
library(knitr)
opts_chunk$set(progress=FALSE, verbose=FALSE, warning=FALSE, message=FALSE, width=1500)

library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)
library(plyr)
library(phangorn)
library(ape)


# Read locations of input files from local umberjack_unit_test.config file
CONFIG_FILENAME <- "./tree_stability.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

# final population concatenated from blocks of sequences at different mutation rates
final_popn_tree_filename <- config[config$key=="FINAL_POPN_TREE_FILENAME",]$val
orig_tree_filename <- config[config$key=="ORIG_TREE_FILENAME",]$val
outdir <- config[config$key=="OUTDIR",]$val

#' Original Coalescent Tree
#' ==========================
#' 

orig_tree <- read.tree(orig_tree_filename)
summary(orig_tree)

#' Depth = `r max(node.depth.edgelength(orig_tree))`
#' 
#' Deepest Tip = `r orig_tree$tip.label[which.max(node.depth.edgelength(orig_tree))]`
#' 
#' Total Tree Branch Length = `r sum(orig_tree$edge.length)`
#' 


#' Final Population Tree  (Constrained By Original Coalescent Topology)
#' ==========================
#' 
final_popn_tree <- read.tree(final_popn_tree_filename)
summary(final_popn_tree)

#' Depth = `r max(node.depth.edgelength(final_popn_tree))`
#' 
#' Deepest Tip = `r final_popn_tree$tip.label[which.max(node.depth.edgelength(final_popn_tree))]`
#' 
#' Total Tree Branch Length = `r sum(final_popn_tree$edge.length)`
#' 


#' Difference Between Original and Final Population Tree
#' ==================================================
#' 
#' 
treedist(orig_tree, final_popn_tree)

#+ fig.width=12, fig.height=15
par(mfrow=c(1,2))
plot(orig_tree, cex=0.6, main="Original Coalescent Tree")
plot(final_popn_tree, cex=0.6, main="Final Population Tree")



#mat <- all.equal.phylo(final_popn_tree, orig_tree, use.edge.length=FALSE, use.tip.label=TRUE, index.return=TRUE)
# 
# 
# # length is 1 as expected
# indelible_tree_filename <- '/home/thuy/gitrepo/SlidingWindow/test/simulations/data/umberjack_unittest/1/tree.nwk'
# indelible_tree <- read.tree(indelible_tree_filename)
# node.depth.edgelength(indelible_tree)
# which.max(node.depth.edgelength(indelible_tree))
# max(node.depth.edgelength(indelible_tree))
# sum(indelible_tree$edge.length)
# 
# fasttree_scratch_filename <- "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/umberjack_unittest/mixed/umberjack_unittest.mixed.fromscratch.nwk"
# fasttree_scratch <- read.tree(fasttree_scratch_filename)
# node.depth.edgelength(fasttree_scratch)
# which.max(node.depth.edgelength(fasttree_scratch))
# max(node.depth.edgelength(fasttree_scratch))
# sum(fasttree_scratch$edge.length)
# 
# 
# # how different is from scratch tree to true toplogy?
# # because fasttree has no idea which is the root, we must arbitrarilyi set a root so that they are the same between both trees
# fasttree_scratch_rooted <- root(fasttree_scratch, outgroup="otu1", resolve.root=TRUE)
# fasttree_scratch_rooted <- multi2di(fasttree_scratch_rooted)
# plot(fasttree_scratch)
# plot(fasttree_scratch_rooted)
# 
# window_tree_rooted <- root(orig_tree, outgroup="otu1", resolve.root=TRUE)
# window_tree_rooted <- multi2di(window_tree_rooted)
# plot(window_tree_rooted, cex=0.5)
# 
# # quite a bit of difference = 38, but this may be because we had to expand polytomies.  Even then, when there were polytomies intact,
# # we got a warning message but the symetric difference was similar:
# # symmetric.difference   branch.score.difference           path.difference quadratic.path.difference 
# # 39.00000                  15.84514                  94.07975                1021.29128 
# # Warning message:
# #   In treedist(window_tree_rooted, fasttree_scratch_rooted) :
# #   Trees are not binary!
# treedist(window_tree_rooted, fasttree_scratch_rooted)
# 
# sum(!window_tree_rooted$tip.label %in% fasttree_scratch_rooted$tip.label )
