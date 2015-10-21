# Checks tree stability from sequence subsampling.
# Do we get the same tree as the full population tree minus unsampled reads when we slide windows across the genome?
library(knitr)
opts_chunk$set(progress=TRUE, verbose=TRUE, warning=FALSE, message=FALSE, width=1500)

library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)
library(plyr)
library(phangorn)
library(ape)
#library(phytools)  # to plot internal node numbers

# Read locations of input files from local umberjack_unit_test.config file
CONFIG_FILENAME <- "./tree_stability.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

# final population concatenated from blocks of sequences at different mutation rates
final_popn_tree_filename <- config[config$key=="FINAL_POPN_TREE_FILENAME",]$val
final_popn_free_tree_filename <- config[config$key=="FINAL_POPN_FREE_TREE_FILENAME",]$val
orig_tree_filename <- config[config$key=="ORIG_TREE_FILENAME",]$val
indelible_tree_filename <- config[config$key=="INDELIBLE_TREE_FILENAME",]$val
outdir <- config[config$key=="OUTDIR",]$val


# Stolen from RF.dist():  reorder tip labels in 1 tree in same order as the other.  
reorder_tips <- function(tree1, tree2) {
  
  ind <- match(tree1$tip.label, tree2$tip.label)
  if (any(is.na(ind)) | length(tree1$tip.label) != length(tree2$tip.label)) 
    stop("trees have different labels")
  tree2$tip.label <- tree2$tip.label[ind]
  ind2 <- match(1:length(ind), tree2$edge[, 2])
  tree2$edge[ind2, 2] <- order(ind)
  return (tree2)
}

# Check if subtrees in tree1 exist in tree2
cmp_subtrees <- function(tree1, tree2) {
  if (is.rooted(tree1)) {
    tree1 <- unroot(tree1)
  }
  
  if (is.rooted(tree2)) {
    tree2 <- unroot(tree2)
  }
    
  sub_t1s <- subtrees(tree1)
  sub_t2s <- subtrees(tree2)
  
  # sort subtrees by number of tips, ascending order
  size_sub_t1s <- sapply(sub_t1s, function(t){length(t$tip)})
  sort_sub_t1s <- sub_t1s[order(size_sub_t1s)]
  
  size_sub_t2s <- sapply(sub_t2s, function(t){length(t$tip)})
  sort_sub_t2s <- sub_t2s[order(size_sub_t2s)]
  
  # Is the subtree in tree1 in tree2?
  not_found <- list()
  for (sub_t1 in sort_sub_t1s) {
    is_found = FALSE
    for (sub_t2 in sort_sub_t2s) {
      if (length(sub_t2$tip) > length(sub_t1$tip)) {  
        break
      } else if (length(sub_t2$tip) == length(sub_t1$tip) &
                   sum (match(sub_t1$tip.label, sub_t2$tip.label, nomatch=0) == 0) == 0) {
        is_found = TRUE
        break          
      }
    }
    
    if (!is_found) {
      not_found <- append(not_found, list(sub_t1))
      
    }
  }
  
  return (not_found)
}

# highlight subtrees in tree1 not found in tree2 as red
plot_tree_diff <- function(tree1, tree2, main) {
  subt1_notin_t2 <- cmp_subtrees(tree1, tree2)
  
  subt1_edges <- unique(unlist(
    sapply(subt1_notin_t2, 
           function(subt1) { which.edge(tree1, subt1$tip.label) })))
  
  tree1_edge_color <- rep("black", length(tree1$edge))
  tree1_edge_color[subt1_edges] <- "red"
  
  subt1_tips <- unique(unlist(sapply(subt1_notin_t2, function(t){t$tip.label})))  
  tree1_tip_color <- ifelse(tree1$tip.label %in% subt1_tips, "red", "black")
  
  plot(tree1, cex=1, main=main, tip.color=tree1_tip_color, edge.color=tree1_edge_color, adj=1)
}


#' Original Coalescent Tree
#' ==========================
#' 

orig_tree <- read.tree(orig_tree_filename)
summary(orig_tree)

#' Depth = `r max(node.depth.edgelength(orig_tree))` = `r max(node.depth.edgelength(orig_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r orig_tree$tip.label[which.max(node.depth.edgelength(orig_tree))]`
#' 
#' Total Tree Branch Length = `r sum(orig_tree$edge.length)`
#' 
#' Depth Summary:
#' 
summary(node.depth.edgelength(orig_tree))


#' Final Population Tree  (Constrained By Original Coalescent Topology)
#' ==========================
#' 
final_popn_tree <- read.tree(final_popn_tree_filename)
summary(final_popn_tree)

#' Depth = `r max(node.depth.edgelength(final_popn_tree))` = `r max(node.depth.edgelength(final_popn_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r final_popn_tree$tip.label[which.max(node.depth.edgelength(final_popn_tree))]`
#' 
#' Total Tree Branch Length = `r sum(final_popn_tree$edge.length)`
#' 
#' Depth Summary:
#' 
summary(node.depth.edgelength(final_popn_tree))


#' Difference Between Original and Final Population Tree
#' ==================================================
#' 
#' 
# to calculate the robinson foulds distance, the trees are unrooted.
treedist(orig_tree, final_popn_tree)


#+ fig.width=20, fig.height=15, dpi=600
par(mfrow=c(1,2), mar=c(0, 0, 4.1, 0))
plot_tree_diff(orig_tree, final_popn_tree, "Original Coalescent Tree")
plot_tree_diff(final_popn_tree, orig_tree, "Final Population Tree")


#' Difference Between Original and Final Population Tree - Remove Duplicate Seq
#' ==================================================
#' 
final_popn_nodup_tree_filename <- config[config$key=="FINAL_POPN_NODUP_TREE_FILENAME",]$val
orig_nodup_tree_filename <- config[config$key=="ORIG_NODUP_TREE_FILENAME",]$val
final_popn_nodup_tree <- read.tree(final_popn_nodup_tree_filename)
orig_nodup_tree <- read.tree(orig_nodup_tree_filename)
treedist(orig_nodup_tree, final_popn_nodup_tree)

#+ fig.width=20, fig.height=15, dpi=600
par(mfrow=c(1,2), mar=c(0, 0, 4.1, 0))
plot_tree_diff(orig_nodup_tree, final_popn_nodup_tree, "Original Coalescent Tree - Duplicate Seq Removed")
plot_tree_diff(final_popn_nodup_tree, orig_nodup_tree, "Final Population Tree - Duplicate Seq Removed")


#' INDELible Tree
#' ================================================
#' 
indelible_tree <- read.tree(indelible_tree_filename)
summary(indelible_tree)

#' Depth = `r max(node.depth.edgelength(indelible_tree))` = `r max(node.depth.edgelength(indelible_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r indelible_tree$tip.label[which.max(node.depth.edgelength(indelible_tree))]`
#' 
#' Total Tree Branch Length = `r sum(indelible_tree$edge.length)`
#' 
#' Depth Summary:
#' 
summary(node.depth.edgelength(indelible_tree))


#' Difference Between Original and INDELible Tree
#' ==================================================
#' 
treedist(orig_tree, indelible_tree)

#+ fig.width=15, fig.height=15
par(mfrow=c(1,2))
plot_tree_diff(orig_tree, indelible_tree, "Original Coalescent Tree")
plot_tree_diff(indelible_tree, orig_tree, "INDELible Tree")


#' Final Population Tree  (No Topology Constraints)
#' ==========================
#' 
final_popn_free_tree <- read.tree(final_popn_free_tree_filename)
summary(final_popn_free_tree)

#' Depth = `r max(node.depth.edgelength(final_popn_free_tree))` = `r max(node.depth.edgelength(final_popn_free_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r final_popn_free_tree$tip.label[which.max(node.depth.edgelength(final_popn_free_tree))]`
#' 
#' Total Tree Branch Length = `r sum(final_popn_free_tree$edge.length)`
#' 
#' Depth Summary:
#' 
summary(node.depth.edgelength(final_popn_free_tree))


#' Difference Between Original and Final Population Tree  (No Topology Constraints)
#' ==================================================
#' 
treedist(orig_tree, final_popn_free_tree)

#+ fig.width=20, fig.height=15
par(mfrow=c(1,2))
plot_tree_diff(orig_tree, final_popn_free_tree, "Original Coalescent Tree")
plot_tree_diff(final_popn_free_tree, orig_tree, "Final Population Tree Made Without Topology Constraints")


#' Difference Between Original and Final Population Tree  (No Topology Constraints) - Duplicate Seq Removed
#' ==================================================
#' 
final_popn_free_nodup_tree_filename <- config[config$key=="FINAL_POPN_FREE_NODUP_TREE_FILENAME",]$val
final_popn_free_nodup_tree <- read.tree(final_popn_free_nodup_tree_filename)
treedist(orig_nodup_tree, final_popn_free_nodup_tree)

#+ fig.width=20, fig.height=15
par(mfrow=c(1,2))
plot_tree_diff(orig_nodup_tree, final_popn_free_nodup_tree, "Original Coalescent Tree - Dup Seq Removed")
plot_tree_diff(final_popn_free_nodup_tree, orig_nodup_tree, "Final Population Tree Made Without Topology Constraints - Dup Seq Removed")



#' Difference Between Final Population Tree Constrained by Original Coalescent Topology and Final Population Tree With No Topology Constraints
#' ==================================================
#' 
treedist(final_popn_tree, final_popn_free_tree)

#+ fig.width=15, fig.height=15
par(mfrow=c(1,2))
plot_tree_diff(final_popn_tree, final_popn_free_tree, main="Final Population Tree Constrained by Original Coalescent Topology")
plot_tree_diff(final_popn_free_tree, final_popn_tree, "Final Population Tree Made Without Topology Constraints")


#' Subsampled Final Population Tree, No Topology Constraints, Dup Seq Removed
#' ======================================================
#' 
final_popn_subsample_free_nodup_tree_filename <- config[config$key=="FINAL_SUBSAMPLE_POPN_FREE_NODUP_TREE_FILENAME",]$val
final_popn_subsample_free_nodup_tree <- read.tree(final_popn_subsample_free_nodup_tree_filename)

#' Depth = `r max(node.depth.edgelength(final_popn_subsample_free_nodup_tree))` = `r max(node.depth.edgelength(final_popn_subsample_free_nodup_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r final_popn_subsample_free_nodup_tree$tip.label[which.max(node.depth.edgelength(final_popn_subsample_free_nodup_tree))]`
#' 
#' Total Tree Branch Length = `r sum(final_popn_subsample_free_nodup_tree$edge.length)`
#' 
#' Depth Summary:
#' 
summary(node.depth.edgelength(final_popn_subsample_free_nodup_tree))


#' Reconstructed Tree From Subsampled Final Population Fasta, No Topology Constraints, Dup Seq Removed
#' ======================================================
#' 
recon_final_popn_subsample_nodup_tree_filename <- config[config$key=="RECON_FINAL_SUBSAMPLE_POPN_NODUP_TREE_FILENAME",]$val
recon_final_popn_subsample_nodup_tree <- read.tree(recon_final_popn_subsample_nodup_tree_filename)

#' Depth = `r max(node.depth.edgelength(recon_final_popn_subsample_nodup_tree))` = `r max(node.depth.edgelength(recon_final_popn_subsample_nodup_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r recon_final_popn_subsample_nodup_tree$tip.label[which.max(node.depth.edgelength(recon_final_popn_subsample_nodup_tree))]`
#' 
#' Total Tree Branch Length = `r sum(recon_final_popn_subsample_nodup_tree$edge.length)`
#' 
#' Depth Summary:
#' 
summary(node.depth.edgelength(recon_final_popn_subsample_nodup_tree))


#' Difference Between Subsampled Final Population Tree and Reconstructed Tree from Subsampled Final Population Seq, Remove Dup Seq
#' ==================================================
#' 

treedist(final_popn_subsample_free_nodup_tree, recon_final_popn_subsample_nodup_tree)

#+ fig.width=20, fig.height=15
par(mfrow=c(1,2))
plot_tree_diff(final_popn_subsample_free_nodup_tree, recon_final_popn_subsample_nodup_tree,main="Subsampled Final Population Tree, No Topology Constraint, Dup Removed")
plot_tree_diff(recon_final_popn_subsample_nodup_tree, final_popn_subsample_free_nodup_tree, main="Reconstructed Tree From Subsampled Final Population Seq, Dup Removed")


#' Resampled Final Population Tree  (Sampling with replacement), No Topology Constraint, Dup Seq Removed
#' ======================================================
#' 
final_popn_resample_free_nodup_tree_filename <- config[config$key=="FINAL_RESAMPLE_POPN_FREE_NODUP_TREE_FILENAME",]$val
final_popn_resample_free_nodup_tree <- read.tree(final_popn_resample_free_nodup_tree_filename)

#' Depth = `r max(node.depth.edgelength(final_popn_resample_free_nodup_tree))` = `r max(node.depth.edgelength(final_popn_resample_free_nodup_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r final_popn_resample_free_nodup_tree$tip.label[which.max(node.depth.edgelength(final_popn_resample_free_nodup_tree))]`
#' 
#' Total Tree Branch Length = `r sum(final_popn_resample_free_nodup_tree$edge.length)`
#' 
#' Depth Summary:
#'
summary(node.depth.edgelength(final_popn_resample_free_nodup_tree))


#' Reconstructed Tree From Resampled Final Population Seq, Dup Seq Removed
#' ===========================================================
#' 
recon_final_popn_reample_nodup_tree_filename <- config[config$key=="RECON_FINAL_RESAMPLE_POPN_NODUP_TREE_FILENAME",]$val
recon_final_popn_resample_nodup_tree <- read.tree(recon_final_popn_reample_nodup_tree_filename)

#' Depth = `r max(node.depth.edgelength(recon_final_popn_resample_nodup_tree))` = `r max(node.depth.edgelength(recon_final_popn_resample_nodup_tree))/1e-4/365`  years of HIV evolution
#' 
#' Deepest Tip = `r recon_final_popn_resample_nodup_tree$tip.label[which.max(node.depth.edgelength(recon_final_popn_resample_nodup_tree))]`
#' 
#' Total Tree Branch Length = `r sum(recon_final_popn_resample_nodup_tree$edge.length)`
#' 
#' Depth Summary:
#'
summary(node.depth.edgelength(recon_final_popn_resample_nodup_tree))



#' Difference Between Resampled Final Population Tree and Reconstructed Tree from Resampled Final Population Seq, Dup Removed
#' ==================================================
#' 


treedist(final_popn_resample_free_nodup_tree, recon_final_popn_resample_nodup_tree)

#+ fig.width=20, fig.height=15
par(mfrow=c(1,2))
plot_tree_diff(final_popn_resample_free_nodup_tree, recon_final_popn_resample_nodup_tree,main="Resampled Final Population Tree, No Topology Constraint, Dup Removed")
plot_tree_diff(recon_final_popn_resample_nodup_tree, final_popn_resample_free_nodup_tree, main="Reconstructed Tree From Resampled Final Population Seq, Dup Removed")
