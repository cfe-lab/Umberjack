#! /usr/bin/Rscript --vanilla

## Quick script to take in a tree with polytomies and use multi2di to
## resolve them.

library(ape);
library(optparse);

option.list <-
  list(make_option("--newick", action="store_true", default=FALSE,
                   help="Parse tree file as Newick [default NEXUS]"),
       make_option("--outnewick", action="store_true", default=FALSE,
                   help="Write output tree file as Newick [default NEXUS]"))

parser <-
  OptionParser(usage="usage: %prog [--newick] [--outnewick] <trees to resolve> <resolved trees>",
               option_list=option.list,
               add_help_option=TRUE)

opts <- parse_args(parser, args=commandArgs(TRUE),
                   positional_arguments=TRUE)

## print(opts)

trees <- NULL;
if (opts$options$newick)
  {
    trees <- read.tree(opts$args[1], keep.multi=TRUE);
  } else {
    trees <- read.nexus(opts$args[1]);
    ## Make sure we have a multiPhylo object coming out of here.
    if (class(trees) == "phylo")
      {
        trees <- c(trees);
      }
  }

resolved.trees <- NULL;
for (i in 1:length(trees))
  {
    tree <- trees[[i]];
    if (i == 1)
      {
        resolved.trees <- c(multi2di(tree));
      } else {
        resolved.trees <- c(resolved.trees, c(multi2di(tree)));
      }
    names(resolved.trees)[i] <- names(trees)[i];
  }

# write.tree takes either a phylo or multiPhylo object; write.nexus
# takes either a phylo or a list of phylos, but it will still work
# with a multiPhylo.
if (opts$options$outnewick)
  {
    write.tree(resolved.trees, file=opts$args[2], tree.names=TRUE);
  } else {
    write.nexus(resolved.trees, file=opts$args[2]);
  }

