#+ setup, include=FALSE
library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, width=1200)


# From data from all windows, aggregates by averaging over windows
library(ggplot2)
library(reshape2)
library(epiR)
library(plyr)
library(stats)
library(psych)
library(GGally)
library(MASS)
NUC_PER_CODON <- 3


# Read locations of input files from local sliding_window_tree_unit_test.config file
CONFIG_FILENAME <- "./aggreg_window.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))
COLLATE_DNDS_FILENAME <- config[config$key=="COLLATE_DNDS_FILENAME",]$val
EXPECTED_DNDS_FILENAME <-  config[config$key=="EXPECTED_DNDS_FILENAME",]$val
EXPECTED_DNDS_START_NUC_POS <-  as.numeric(config[config$key=="EXPECTED_DNDS_START_NUC_POS",]$val)
EXPECTED_DNDS_END_NUC_POS <-  as.numeric(config[config$key=="EXPECTED_DNDS_END_NUC_POS",]$val)


FULL_POPN_CONSERVE_CSV <- config[config$key=="FULL_POPN_CONSERVE_CSV",]$val
ORIG_CONSERVE_CSV <- config[config$key=="ORIG_CONSERVE_CSV",]$val
ALN_CONSERVE_CSV <- config[config$key=="ALN_CONSERVE_CSV",]$val

SMOOTH_DIST <-  as.numeric(config[config$key=="SMOOTH_DIST",]$val)

#+ results='asis'
kable(config, format="html", caption="config")

#+
# Cols:  Window_Start,Window_End,Reads,CodonSite,CodonDepth,Conserve,Entropy,N,S,EN,ES,dN,dS,dN_minus_dS
collate_dnds <- read.table(COLLATE_DNDS_FILENAME, header=TRUE, na.strings="None", comment.char = "#", sep=",")
collate_dnds$Subst <- collate_dnds$N + collate_dnds$S
collate_dnds$Diverge <- 1 - collate_dnds$Conserve
colnames(collate_dnds)[grep("^Codons$", colnames(collate_dnds), perl=TRUE)] <- "CodonDepth"
collate_dnds$dNdS <- collate_dnds$dN/collate_dnds$dS
collate_dnds$dNdS[collate_dnds$dS==0] <- NA
collate_dnds$ErrFractn <- collate_dnds$Err/(collate_dnds$Reads)
collate_dnds$AmbigCodonFractn <- 1 - collate_dnds$CodonDepth/collate_dnds$Reads
collate_dnds$Err_N_Fractn <- collate_dnds$Err_N/collate_dnds$CodonDepth
collate_dnds$Err_S_Fractn <- collate_dnds$Err_S/collate_dnds$CodonDepth 
collate_dnds$Pad_Fractn <- collate_dnds$Pad/(collate_dnds$Reads)

# Average across all codon sites in a window
per_window_ave <- ddply(.data=collate_dnds, .variables="Window_Start", 
                        .fun=function(x) {                            
                          data.frame(Window_Diverge=mean(x$Diverge, na.rm=TRUE),
                                     Window_Conserve=mean(x$Conserve, na.rm=TRUE),
                                     Window_Entropy=mean(x$Entropy, na.rm=TRUE),
                                     Window_Subst=mean(x$Subst, na.rm=TRUE),
                                     Window_CodonDepth=mean(x$CodonDepth, na.rm=TRUE),
                                     Window_ErrFractn=mean(x$ErrFractn, na.rm=TRUE),
                                     Window_AmbigCodonFractn=mean(x$AmbigCodonFractn, na.rm=TRUE),
                                     Window_Err_N_Fractn=mean(x$Err_N_Fractn, na.rm=TRUE),
                                     Window_Err_S_Fractn=mean(x$Err_S_Fractn, na.rm=TRUE)
                                     )
                        })
collate_dnds <- merge(x=collate_dnds, y=per_window_ave, by="Window_Start", all=TRUE, sort=TRUE)

dim(collate_dnds)
head(collate_dnds)

#' **Summary Per-Window-Codon Stats**
#' 
summary(collate_dnds)


#+
# Count windows per Codon Site.  Assume that each window represented once.  Assume each codon site represented once per window.
actual_dnds <- as.data.frame(with(collate_dnds, table(CodonSite)))
colnames(actual_dnds)[grep("Freq", colnames(actual_dnds))] <- "Windows"
actual_dnds$CodonSite <- as.numeric(actual_dnds$CodonSite)  # table converts CodonSite to a factor

# Average over windows
per_codon_ave <- ddply(.data=collate_dnds, .variables="CodonSite", 
                          .fun=function(x) {                            
                            hi_syn_x <- x[x$S >=1, ]
                            hi_subst_x <- x[x$N >=1 & x$S >=1, ]
                            
                            data.frame(
                             aveDnDs=mean(x$dNdS, na.rm=TRUE),
                             CodonDepth=mean(x$CodonDepth, na.rm=TRUE),
                             Conserve=mean(x$Conserve, na.rm=TRUE),
                             Diverge=mean(1-x$Conserve, na.rm=TRUE),
                             Entropy=mean(x$Entropy, na.rm=TRUE),
                             Reads=mean(x$Reads, na.rm=TRUE), 
                             NonSyn=mean(x$N, na.rm=TRUE), 
                             Syn=mean(x$S, na.rm=TRUE), 
                             Subst=mean(x$Subst, na.rm=TRUE), 
                             dN=mean(x$dN, na.rm=TRUE), 
                             dS=mean(x$dS, na.rm=TRUE), 
                             dN_minus_dS=mean(x$dN_minus_dS, na.rm=TRUE),
                             dN_minus_dSWeightByReadsNoLowSubst=weighted.mean(hi_subst_x$dN_minus_dS, hi_subst_x$Reads, na.rm=TRUE),
                             dNdSWeightBySubst=weighted.mean(x$dNdS, x$Subst, na.rm=TRUE),
                             dNdSWeightByReads=weighted.mean(x$dNdS, x$Reads, na.rm=TRUE),
                             dNdSWeightByReadsNoLowSyn=weighted.mean(hi_syn_x$dNdS, hi_syn_x$Reads, na.rm=TRUE),
                             dNdSWeightByReadsNoLowSubst=weighted.mean(hi_subst_x$dNdS, hi_subst_x$Reads, na.rm=TRUE),
                             dN_minus_dSWeightByReads=weighted.mean(x$dN_minus_dS, x$Reads, na.rm=TRUE))
                        })
actual_dnds <- merge(x=actual_dnds, y=per_codon_ave, by="CodonSite", all=TRUE, sort=TRUE)

#+
# Get average smoothed over multiple sites
MIN_COLLATE_SITE <- min(collate_dnds$CodonSite)
MAX_COLLATE_SITE <- max(collate_dnds$CodonSite)

actual_dnds$multisiteAvedNdS <- sapply(actual_dnds$CodonSite, function(site) {
  mean(collate_dnds[collate_dnds$CodonSite >= (site - SMOOTH_DIST) & collate_dnds$CodonSite <= (site + SMOOTH_DIST), "dNdS"], 
       na.rm=TRUE)
})

actual_dnds$multisiteAvedNdSWeightBySubst <- sapply(actual_dnds$CodonSite, function(site) {  
  set <- collate_dnds[collate_dnds$CodonSite >= (site - SMOOTH_DIST) & 
                              collate_dnds$CodonSite <= (site + SMOOTH_DIST),]
  weighted.mean(set$dNdS, set$Subst, na.rm=TRUE)
})

# Here, we don't calculate the dN/dS ratio until the very end.
# We are naturally weighting the sites by their number of substitutions
actual_dnds$multisitedNdSWeightBySubst <- sapply(actual_dnds$CodonSite, function(site) {
  set <- collate_dnds[collate_dnds$CodonSite >= (site - SMOOTH_DIST) & 
                        collate_dnds$CodonSite <= (site + SMOOTH_DIST),]
  total_nonsyn <- sum(set$N, na.rm=TRUE)
  total_syn <- sum(set$S, na.rm=TRUE)
  total_exp_nonsyn <- sum(set$EN, na.rm=TRUE)
  total_exp_syn <- sum(set$ES, na.rm=TRUE)
  multisite_dnds <- (total_nonsyn * total_exp_syn) / (total_syn * total_exp_nonsyn )
})


#' **Summary Per-Codon Stats (aggregated over windows)**
#' 
summary(actual_dnds)



#'
#'  Read in Expected dN/dS
#'  ==============================
#'  

# Cols: Observed S Changes  Observed NS Changes  E[S Sites]	E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS	dn/ds
# Parse the expected dnds filename for it start and end nucleotide positions (1-based)
expected_dnds <- read.table(EXPECTED_DNDS_FILENAME, header=TRUE, sep="\t")  # Site  Interval	Scaling_factor	Rate_class	Omega
expected_dnds_start_codon <- (EXPECTED_DNDS_START_NUC_POS %/% 3) + 1
expected_dnds$CodonSite <- as.numeric(rownames(expected_dnds)) + expected_dnds_start_codon -1
expected_dnds$Omega <- expected_dnds$dN/expected_dnds$dS
expected_dnds$Omega[expected_dnds$dS == 0] <- NA
colnames(expected_dnds)[grep("Observed.S.Changes", colnames(expected_dnds))] <- "S"
colnames(expected_dnds)[grep("Observed.NS.Changes", colnames(expected_dnds))] <- "N"
expected_dnds$Subst <- expected_dnds$S + expected_dnds$N


# check consistency
if (!all(expected_dnds$CodonSite == actual_dnds$CodonSite)) {
  stop("Expected dN/dS Sites don't line up with Actual dN/dS Sites")
}

# TODO:  don't need the for loops
MIN_SITE <- min(expected_dnds$CodonSite)
MAX_SITE <- max(expected_dnds$CodonSite)
expected_dnds$MultisiteAveDnDs <- apply(expected_dnds, 1, function(row) {
  site <- as.numeric(row["CodonSite"])
  total_dnds <- 0
  total_sites <- 0
  for (adjsite in max(MIN_SITE, site-SMOOTH_DIST):min(MAX_SITE, site+SMOOTH_DIST)) {
    if (!is.na(expected_dnds[expected_dnds$CodonSite==adjsite, "Omega"])) {
      total_dnds <- total_dnds + expected_dnds[expected_dnds$CodonSite==adjsite, "Omega"]
      total_sites <- total_sites +1
    }
  }
  multisite_ave <- total_dnds / total_sites
  return (multisite_ave)
})

#' **Summary Expected dN/dS Stats**
#' 
summary(expected_dnds)



#' **How is our coverage?**
#' -------------------------------------

# Plots the window stats for each window at each codon site
plot_window_nums <- function(collate_dnds, colname, descname) {
  fig <- ggplot(collate_dnds, aes(x=CodonSite)) + 
    geom_point( aes_string(y=colname), shape=1, alpha=0.5) +  
    geom_smooth( aes_string(y=colname)) + 
    ylab(paste0(descname, "\n")) + 
    xlab("\nCodon Site") + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0(descname, " From Each Window"))
  print(fig)
}


plot_window_nums(collate_dnds, "Reads", "Max Read Depth")
plot_window_nums(collate_dnds, "CodonDepth", "Unambiguous Codon Count")
plot_window_nums(collate_dnds, "Diverge", "Divergence")
plot_window_nums(collate_dnds, "Entropy", "Entropy")
plot_window_nums(collate_dnds, "Subst", "Substitutions")
plot_window_nums(collate_dnds, "N", "Nonsynonymous Substitutions")
plot_window_nums(collate_dnds, "S", "Synonymous Substitutions")


#' How much variation is there from window to window for the same codon site?
#' -------------------------------------------------------------------------------
#' 

# Plots the coefficient of variance across windows at each codon site
# Assumes that collate_dnds has data
plot_coeff_var <- function(collate_dnds, colname, descname) {
  coeff_var <- ddply(.data=collate_dnds, 
                     .variables="CodonSite",                          
                     .fun=function(x) {
                           ave <- mean(x[, colname], na.rm=TRUE)
                           stddev <- sd(x[, colname], na.rm=TRUE)
                           cv <- stddev/ave
                           return (data.frame(ave=ave, stddev=stddev, cv=cv))
                         })
  
  coeff_var <- coeff_var[order(coeff_var$CodonSite), ]
  
  # Check that we have the correct dimensions
  if(nrow(coeff_var) != length(unique(collate_dnds$CodonSite))) {
    stop("Invalid dimensions.  Number of codon sites for coeff_var != number of codon sites for collate_dnds")
  }
  
  print(summary(coeff_var))

  fig <- ggplot(coeff_var, aes(x=CodonSite, y=cv)) + 
    geom_point(shape=1, alpha=0.5) +  
    geom_smooth() + 
    ylab(paste0("Coeff of Var of ", descname, "\n")) + 
    xlab("\nCodon Site") + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Coefficient of Variance Across Windows for ", descname))
  print(fig)
}
  
plot_coeff_var(collate_dnds, "Conserve", "Conservation")
plot_coeff_var(collate_dnds, "Reads", "Max Reads")
plot_coeff_var(collate_dnds, "Diverge", "Divergence")
plot_coeff_var(collate_dnds, "Entropy", "Entropy")
plot_coeff_var(collate_dnds, "N", "Nonsynonymous Substitutions")
plot_coeff_var(collate_dnds, "S", "Synonymous Substitutions")
plot_coeff_var(collate_dnds, "Subst", "Substitutions")
plot_coeff_var(collate_dnds, "dNdS", "dN/dS")
plot_coeff_var(collate_dnds, "dN", "dN")
plot_coeff_var(collate_dnds, "dS", "dS")
plot_coeff_var(collate_dnds, "dN_minus_dS", "Scaled dN-dS")


#' What is the Relationship Between Sequence Diversity And Full Population Phylogeny Substitution Count?
#' ------------------------------------------
#' 

# Takes a *.conserve.csv file with nucleotide stats for unsliced sequences
# Expects cols: "NucSite", "Conserve", "Entropy", "NucDepth", "CodonDepth"  
# Returns a dataframe with elements:
# - CodonDat = dataframe with values average per-codon values for Conserve, Diverge, Entropy, CodonDepth
# - NucStats = matrix summary of per-nucleotide stats
fill_nowindow_codon_conserve_dat <- function(nuc_conserve_csv) {
  # Per-nucleotide conservation, Divergence, Entropy
  conserve_nuc <- read.table(nuc_conserve_csv, header=TRUE, sep=",")
  conserve_nuc$Diverge <- 1 - conserve_nuc$Conserve
  
  NUM_NUC_SITES <- max(conserve_nuc$NucSite)
  NUM_CODON_SITES <- ceiling(NUM_NUC_SITES / 3)
  
  # Get per-codon average stats for Divergence, Entropy
  conserve_codon <- adply(.data=c(1:NUM_CODON_SITES), .margins=1,
                           .fun=function(codonsite) {
                             start_nuc_site <- (codonsite * 3) -2  # 1-based nucleotide sites
                             end_nuc_site <- codonsite * 3
                             codon_set <- conserve_nuc[conserve_nuc$NucSite >= start_nuc_site & conserve_nuc$NucSite <= end_nuc_site,]
                             data.frame(  
                               CodonDepth=conserve_nuc[start_nuc_site, "CodonDepth"],
                               Conserve=mean(codon_set$Conserve),
                               Diverge=mean(codon_set$Diverge),
                               Entropy=mean(codon_set$Entropy))
                          })
  colnames(conserve_codon)[1] <- "CodonSite"  # rename X1 column to CodonSite
  conserve_codon$CodonSite <- as.numeric(conserve_codon$CodonSite)
  returnDat <- list(CodonDat=conserve_codon,
                          NucStat=summary(conserve_nuc))
  return (returnDat)
}



# Scatterplots Phylogeny Substitution Count Vs Sequence Diversity
# INPUT:
# - codon_dat:  dataframe with cols: CodonSite,CodonDepth,Entropy,N,S,Subst,Diverge
# - diversity_title:  x-axes title prefix
# - subst_title:  y-axes title prefix
# - diversity_measure_suffix:  diversity measure suffix
# - subst_measure_suffix:  subsitution measure suffix
# - codon_diversity:  dataframe with cols: CodonSite,CodonDepth,Entropy,Diverge.  Only used if codon_dat is NULL
# - codon_diversity:  dataframe with cols: CodonSite,N,S, Subst.  Only used if codon_dat is NULL
plot_subst_vs_diversity <- function(codon_dat, diversity_title, subst_title, diversity_measure_suffix, subst_measure_suffix, 
                                    codon_diversity=NULL, codon_subst=NULL) {
  if (is.null(codon_dat) & (is.null(codon_diversity)  | is.null(codon_subst))) {
    stop("Either define codon_dat or both codon_diversity and codon_subst")
  } else if (is.null(codon_dat))  {
    codon_dat <- merge(x=subset(codon_subst, select=c(CodonSite, S, N, Subst)),
                       y=subset(codon_diversity, select=c(CodonSite, CodonDepth, Diverge, Entropy)),
                       by="CodonSite", all=TRUE)
  }
  phylo_subst_melt <- reshape2:::melt(data=subset(codon_dat, select=c(CodonSite, S, N, Subst)),
                                             id.vars="CodonSite",
                                             measure.vars=c("S", "N", "Subst"),
                                             variable.name="PhyloMeasure", 
                                             value.name="PhyloVal")
  phylo_subst_melt$PhyloMeasure <- paste0(phylo_subst_melt$PhyloMeasure, subst_measure_suffix)
  
  leaf_diversity_melt <- reshape2:::melt(data=subset(codon_dat, select=c(CodonSite, CodonDepth, Diverge, Entropy)),
                                                id.vars="CodonSite",
                                                measure.vars=c("CodonDepth", "Diverge", "Entropy"),
                                                variable.name="LeafMeasure", 
                                                value.name="LeafVal")
  leaf_diversity_melt$LeafMeasure <- paste0(leaf_diversity_melt$LeafMeasure, diversity_measure_suffix)
  
  combo <- merge(x=phylo_subst_melt, y=leaf_diversity_melt, by="CodonSite", all=TRUE)
  
  fig <- ggplot(combo, aes(x=LeafVal, y=PhyloVal)) + 
    facet_grid(PhyloMeasure~LeafMeasure, scales="free_x") + 
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm") + 
    xlab(paste0("\n", diversity_title, " Sequence Diversity")) + 
    ylab(paste0(subst_title, " SLAC Count\n")) + 
    ggtitle(paste0("Scatterplot ", subst_title, " SLAC Count Vs ", diversity_title, " Sequence Diversity"))
  print(fig)
}

#'
#' Window Sequence Diversity Vs Window Phylogeny Substitution Count
#' =====================================================================
#'  

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(collate_dnds, diversity_title="Window", subst_title="Window", 
                        diversity_measure_suffix=".Win", subst_measure_suffix=".Win")

#' Full Population Sequence Diversity Vs Full Population Phylogeny Substitution Count
#' =====================================================================
#'  
full_popn_conserve_dat <- fill_nowindow_codon_conserve_dat(FULL_POPN_CONSERVE_CSV)
#+ results='asis'
kable(full_popn_conserve_dat$NucStat, format="html", caption="Full Population Nucleotide Stats")
kable(summary(full_popn_conserve_dat$CodonDat), format="html", caption="Full Population Codon Stats")

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(codon_dat=NULL, diversity_title="Full Population", subst_title="Full Population", 
                        diversity_measure_suffix=".Full", subst_measure_suffix=".Full", 
                        codon_diversity=full_popn_conserve_dat$CodonDat, codon_subst=expected_dnds)


#'
#' Original Read Diversity Vs Full Population Phylogeny Substitution Count
#' =====================================================================
#'  

orig_conserve_dat <- fill_nowindow_codon_conserve_dat(ORIG_CONSERVE_CSV)

#+ results='asis'
kable(orig_conserve_dat$NucStat, format="html", caption="Original Reads Nucleotide Stats")
kable(summary(orig_conserve_dat$CodonDat), format="html", caption="Original Reads Codon Stats")

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(codon_dat=NULL, diversity_title="Original Read", subst_title="Full Population", 
                        diversity_measure_suffix=".OrigRead", subst_measure_suffix=".Full", 
                        codon_diversity=orig_conserve_dat$CodonDat, codon_subst=expected_dnds)

#'
#'  Aligned Read Diversity Vs Full Population Phylogeny Substitution Count
#'  =====================================================================
#'  
aln_conserve_dat <- fill_nowindow_codon_conserve_dat(ALN_CONSERVE_CSV)

#+ results='asis'
kable(aln_conserve_dat$NucStat, format="html", caption="Aligned Reads Nucleotide Stats")
kable(summary(aln_conserve_dat$CodonDat), format="html", caption="Aligned Reads Codon Stats")

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(codon_dat=NULL, diversity_title="Aligned Read", subst_title="Full Population", 
                        diversity_measure_suffix="Aligned", subst_measure_suffix=".Full", 
                        codon_diversity=aln_conserve_dat$CodonDat, codon_subst=expected_dnds)

#'
#' What is the Effect of Window-ing?
#' ------------------------------------------
#' 
# Scatterplots Window Sequence Diversity Vs Non Window Sequence Diversity
# Scaterplots Window SLAC Counts Vs Non Window Sequence Diversity
# Scaterplots Window SLAC Counts Vs Non Window SLAC Counts
# I.e.  Compares effect of Window-ing on Diversity.  Compares effect of Window-ing on Substitutions.
# INPUT:
# - nowindow_leaf_diversity:  dataframe with cols:  CodonSite, CodonDepth, Diverge, Entropy.  
#     Values are calculated per-nucleotide and averaged per-codon.  No windows used in calculations.
# - nowindow_phylo_subst:  dataframe with cols: N, S, Subst
#     Values are calculated per codon against full population.  No windows used in calculations.  No reads used in calculations.
# - win_codon_dat:  dataframe with cols: CodonSite,CodonDepth,Entropy,N,S,Subst,Diverge 
# - nowindow_title:  nowindow axes title
# - window_title:  window_title axes title
# - nowindow_measure_suffix:  nowindow measure suffix
# - window_measure_suffix:  window measure suffix
plot_win_subst_div_vs_nowin_div <- function(nowindow_leaf_diversity, nowindow_phylo_subst, win_codon_dat,
                                     nowindow_title, window_title, nowindow_measure_suffix, window_measure_suffix) {
  
  nowindow_leaf_diversity_melt <- reshape2:::melt(data=subset(nowindow_leaf_diversity, select=c(CodonSite, CodonDepth, Diverge, Entropy)),
                                                  id.vars="CodonSite",
                                                  measure.vars=c("CodonDepth", "Diverge", "Entropy"),
                                                  variable.name="LeafMeasure.NoWin", 
                                                  value.name="LeafVal.NoWin")
  nowindow_leaf_diversity_melt$LeafMeasure.NoWin <- paste0(nowindow_leaf_diversity_melt$LeafMeasure.NoWin, nowindow_measure_suffix)
  
  nowindow_phylo_subst_melt <- reshape2:::melt(data=subset(nowindow_phylo_subst, select=c(CodonSite, S, N, Subst)),
                                               id.vars="CodonSite",
                                               measure.vars=c("S", "N", "Subst"),
                                               variable.name="PhyloMeasure.NoWin", 
                                               value.name="PhyloVal.NoWin")
  nowindow_phylo_subst_melt$PhyloMeasure.NoWin <- paste0(nowindow_phylo_subst_melt$PhyloMeasure.NoWin, nowindow_measure_suffix)
  
  window_phylo_subst_melt <- reshape2:::melt(data=subset(win_codon_dat, select=c(CodonSite, S, N, Subst)),
                                             id.vars="CodonSite",
                                             measure.vars=c("S", "N", "Subst"),
                                             variable.name="PhyloMeasure.Win", 
                                             value.name="PhyloVal.Win")
  window_phylo_subst_melt$PhyloMeasure.Win <- paste0(window_phylo_subst_melt$PhyloMeasure.Win, window_measure_suffix)
  
  window_leaf_diversity_melt <- reshape2:::melt(data=subset(win_codon_dat, select=c(CodonSite, CodonDepth, Diverge, Entropy)),
                                                id.vars="CodonSite",
                                                measure.vars=c("CodonDepth", "Diverge", "Entropy"),
                                                variable.name="LeafMeasure.Win", 
                                                value.name="LeafVal.Win")
  window_leaf_diversity_melt$LeafMeasure.Win <- paste0(window_leaf_diversity_melt$LeafMeasure.Win, window_measure_suffix)
  
  
#   win_nowindow_leaf_combo <- merge(x=window_leaf_diversity_melt, y=nowindow_leaf_diversity_melt, by="CodonSite", all=TRUE)
#   fig <- ggplot(win_nowindow_leaf_combo, aes(x=LeafVal.NoWin, y=LeafVal.Win)) + 
#     facet_grid(LeafMeasure.Win~LeafMeasure.NoWin, scales="free") + 
#     geom_point(shape=1, alpha=0.1) + 
#     geom_smooth(method="lm") + 
#     xlab(paste0("\n", nowindow_title, " Sequence Diversity")) + 
#     ylab(paste0(window_title, " Sequence Diversity\n")) + 
#     ggtitle(paste0("Scatterplot ", window_title, " Sequence Diversity Vs ", nowindow_title, " Sequence Diversity"))
#   print(fig)
#   
  win_nowindow_leaf_combo <- merge(x=subset(win_codon_dat, select=c(CodonSite, CodonDepth, Diverge, Entropy)), 
                                       y=subset(nowindow_leaf_diversity, select=c(CodonSite, CodonDepth, Diverge, Entropy)), 
                                   by="CodonSite", all.x=TRUE,  all.y=FALSE, 
                                   suffixes=c(window_measure_suffix, nowindow_measure_suffix))
  
  sapply(c("CodonDepth", "Diverge", "Entropy"), function(colname) {
      fig <- ggplot(win_nowindow_leaf_combo, aes_string(x=paste0(colname, nowindow_measure_suffix), 
                                                        y=paste0(colname, window_measure_suffix))) + 
        geom_point(shape=1, alpha=0.1) + 
        geom_smooth(method="lm") + 
        xlab(paste0("\n", nowindow_title, " ", colname)) + 
        ylab(paste0(window_title, " ", colname, "\n")) + 
        ggtitle(paste0("Scatterplot ", window_title, " ", colname, " Vs ", nowindow_title, " ", colname))
      print(fig)
  })
  
  # Scatterplot Window Phylogeny Substitutions Vs  Non-Window Sequence Diversity  
  win_subst_nowindow_leaf_combo <- merge(x=window_phylo_subst_melt, y=nowindow_leaf_diversity_melt, by="CodonSite", all=TRUE)
  fig <- ggplot(win_subst_nowindow_leaf_combo, aes(x=LeafVal.NoWin, y=PhyloVal.Win)) + 
    facet_grid(PhyloMeasure.Win~LeafMeasure.NoWin, scales="free_x") + 
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm") + 
    xlab(paste0("\n", nowindow_title, " Sequence Diversity")) + 
    ylab(paste0(window_title, " SLAC Count\n")) + 
    ggtitle(paste0("Scatterplot ", window_title, " SLAC Count Vs ", nowindow_title, " Sequence Diversity"))
  print(fig)
}

# Scatterplots Window SLAC Counts vs Full Population SLAC Counts
plot_win_subst_vs_full_subst <- function(nowindow_phylo_subst, win_codon_dat,
                                            nowindow_title, window_title, nowindow_measure_suffix, window_measure_suffix) {
  
  # Density Plot of Window-Codon Sites at Each Full Population N, S, Subst Count
  # Scatterplot Window SLAC Vs  Full Population SLAC
  win_nowindow_subst_combo <- merge(x=subset(nowindow_phylo_subst, select=c(CodonSite, S, N, Subst)),
                                    y=subset(win_codon_dat, select=c(CodonSite, S, N, Subst)), 
                                    by="CodonSite", all=TRUE, suffixes=c(nowindow_measure_suffix, window_measure_suffix))
  
  sapply(c("N", "S", "Subst"), function(col) {
    
    # Cumulative Frequency plot Window-Codon Sites at Each Full POpulation N, S, Subst Count
    fig <- ggplot(win_nowindow_subst_combo, aes_string(x=paste0(col, nowindow_measure_suffix))) +
        stat_ecdf() + 
        xlab(paste0("\n", nowindow_title, " " , col)) + 
        ylab("Cumulative Frequency Window-Codon Sites\n") + 
        ggtitle(paste0("Cumulative Frequency Plot Window-Codon Sites by  ", nowindow_title, " ", col))
    print(fig)
      
    
    # Scatterplot Window N, S, Subst Vs Full Population N, S, Subst
    fig <- ggplot(win_nowindow_subst_combo, aes_string(x=paste0(col, nowindow_measure_suffix), 
                                                       y=paste0(col, window_measure_suffix))) + 
      geom_point(shape=1, alpha=0.1) + 
      geom_smooth(method="lm") + 
      geom_abline(slope=1, color="red") + 
      xlab(paste0("\n", nowindow_title, " " , col)) + 
      ylab(paste0(window_title, " ", col, "\n")) + 
      ggtitle(paste0("Scatterplot ", window_title, " ", col, " Vs ", nowindow_title, " ", col))
    print(fig)
  })
  
}

#'
#' Effect of Windowing on Phylogeny Compared to Full Population
#' =================================================
#' 

#+ results='asis'
kable(full_popn_conserve_dat$NucStat, format="html", caption="Full Population Nucleotide Stats")
kable(summary(full_popn_conserve_dat$CodonDat), format="html", caption="Full Population Codon Stats")

#+ fig.width=12, fig.height=12
plot_win_subst_vs_full_subst(nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
                                nowindow_title="Full Population", window_title="Window", 
                                nowindow_measure_suffix=".Full", window_measure_suffix=".Win")

#' **What happens to windows when full population N > mean + 1SD = `r mean(expected_dnds$N) + sd(expected_dnds$N)` ?**
#' 
summary(collate_dnds[collate_dnds$CodonSite %in% expected_dnds[expected_dnds$N > (mean(expected_dnds$N) + sd(expected_dnds$N)), 
                                                               "CodonSite"], ])

#' **What happens to windows when full population S > mean+1SD = `r mean(expected_dnds$S) + sd(expected_dnds$S)` ?**
#' 
summary(collate_dnds[collate_dnds$CodonSite %in% expected_dnds[expected_dnds$S > (mean(expected_dnds$S) + sd(expected_dnds$S)), 
                                                               "CodonSite"], ])

#' **What happens to windows when full population Subst > mean+1SD = `r mean(expected_dnds$Subst) + sd(expected_dnds$Subst)` ?**
#' 
summary(collate_dnds[collate_dnds$CodonSite %in% expected_dnds[expected_dnds$Subst > (mean(expected_dnds$Subst) + sd(expected_dnds$Subst)),
                                                               "CodonSite"], ])


#'
#' Effect on Windowing in Sequencing Diversity and Phylogeny Substitutions Compared to Full Population
#' =================================================
#' 

#+ fig.width=12, fig.height=12
plot_win_subst_div_vs_nowin_div(nowindow_leaf_diversity=full_popn_conserve_dat$CodonDat, 
                         nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
                         nowindow_title="Full Population", window_title="Window", 
                         nowindow_measure_suffix=".Full", window_measure_suffix=".Win")
#'
#' Effect on Windowing in Sequencing Diversity and Phylogeny Substitutions Compared to Original ART reads
#' =================================================
#' 

#+ results='asis'
kable(orig_conserve_dat$NucStat, format="html", caption="Original Reads Nucleotide Stats")
kable(summary(orig_conserve_dat$CodonDat), format="html", caption="Original Reads Codon Stats")

#+ fig.width=12, fig.height=12
plot_win_subst_div_vs_nowin_div(nowindow_leaf_diversity=orig_conserve_dat$CodonDat, 
                         nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
                         nowindow_title="Original Reads", window_title="Window", 
                         nowindow_measure_suffix=".Orig", window_measure_suffix=".Win")


#' Effect on Windowing in Sequencing Diversity and Phylogeny Substitutions Compared to Aligned ART reads
#' =================================================
#' 

#+ results='asis'
kable(aln_conserve_dat$NucStat, format="html", caption="Aligned Reads Nucleotide Stats")
kable(summary(aln_conserve_dat$CodonDat), format="html", caption="Aligned Reads Codon Stats")

#+ fig.width=12, fig.height=12
plot_win_subst_div_vs_nowin_div(nowindow_leaf_diversity=aln_conserve_dat$CodonDat, 
                                nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
                                nowindow_title="Aligned Reads", window_title="Window", 
                                nowindow_measure_suffix=".Aln", window_measure_suffix=".Win")




#'  Compare dN/dS against Expected
#'  -------------------------------------------------------
#'  

#' **Scatterplot actual vs expected dn ds together**
fullDat <- cbind(actual_dnds, 
                 Expected=expected_dnds$Omega, 
                 ExpectedMultisite=expected_dnds$MultisiteAveDnDs, 
                 ExpectedMinus=expected_dnds$Scaled.dN.dS,
                 ExpectedSyn=expected_dnds$S,
                 ExpectedNonSyn=expected_dnds$N)
summary(fullDat)


# Expects that fullDat has been filled
scatterplot_act_exp_dnds <- function(act_colname, act_descname, exp_colname, exp_descname) {
  fig <- ggplot(fullDat , aes_string(x=exp_colname, y=act_colname)) + 
    geom_point() +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0(exp_descname, " Vs ", act_descname))
  print(fig)
}

#'Scatterplot Expected vs Inferred Averaged Across Windows
#'  
#+ fig.width=10, fig.height=10
scatterplot_act_exp_dnds("aveDnDs", "Inferred dN/dS Averaged Across Windows", "Expected", "Expected dN/dS")
scatterplot_act_exp_dnds("dNdSWeightBySubst", 
                         "Inferred dN/dS Averaged Across Windows, Weighted by Substitutions", 
                         "Expected", "Expected dN/dS")
scatterplot_act_exp_dnds("dNdSWeightByReads", "Inferred dN/dS Averaged Across Windows, Weighted by Reads", "Expected", "Expected dN/dS")
scatterplot_act_exp_dnds("dNdSWeightByReadsNoLowSubst", 
                         "Inferred dN/dS Averaged Across Windows With >=1 N, >=1 S, Weighted by Reads", 
                         "Expected", "Expected dN/dS")
scatterplot_act_exp_dnds("dNdSWeightByReadsNoLowSyn", 
                         "Inferred dN/dS Averaged Across Windows With >=1 S, Weighted by Reads", 
                         "Expected", "Expected dN/dS")
scatterplot_act_exp_dnds("dN_minus_dS", "Inferred Scaled dN-dS Averaged Across Windows", "ExpectedMinus", "Expected Scaled dN-dS")
scatterplot_act_exp_dnds("dN_minus_dSWeightByReadsNoLowSubst", 
                         "Inferred Scaled dN-dS Averaged Across Windows With >=1 N, >=1 S, Weighted by Reads", 
                         "ExpectedMinus", "Expected Scaled dN-dS")
scatterplot_act_exp_dnds("multisiteAvedNdS", "Inferred Multisite Ave dN/dS", 
                         "ExpectedMultisite", "Expected Multisite Ave dN/dS")
scatterplot_act_exp_dnds("multisiteAvedNdSWeightBySubst", "Inferred Multisite Ave dN/dS, Sites weighted by Substitutions", 
                         "ExpectedMultisite", "Expected Multisite Ave dN/dS")
scatterplot_act_exp_dnds("multisitedNdSWeightBySubst", "Inferred Multisite dN/dS", 
                         "ExpectedMultisite", "Expected Multisite Ave dN/dS")


fullDatBydatsource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars=c("CodonSite"),
                                                 variable.name="datsource", value.name="perSiteVal")
head(fullDatBydatsource)
tail(fullDatBydatsource)
str(fullDatBydatsource)
summary(fullDatBydatsource)

# Plots the dn/ds across the genome
# Assumes that fullDatBydatsource has been filled in
plot_dnds_across_genome <- function(cols, title, ylabel) {
  # Smoothed 
  fig <- ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% cols,], 
         aes(x=CodonSite, y=perSiteVal, color=datsource, linetype=datsource) ) + 
    geom_smooth() + 
    xlab("\nCodon Site") + 
    ylab(paste0(ylabel, "\n")) + 
    ggtitle(paste0("Smoothed ", title, " by Codon Site")) + 
    theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
          legend.text=element_text(size=24), legend.title=element_blank())
  print(fig)
  
  # Points
  fig <- ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% cols,], 
                aes(x=CodonSite, y=perSiteVal, color=datsource) ) + 
    geom_line(shape=1, alpha=0.5) + 
    xlab("\nCodon Site") + 
    ylab(paste0(ylabel, "\n")) + 
    ggtitle(paste0(title, " by Codon Site")) + 
    theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
          legend.text=element_text(size=24), legend.title=element_blank())
  print(fig)
}


#' **Smoothed Scatterplot of Site dn/ds across the genome**
#' 
#+ fig.width=20
plot_dnds_across_genome(cols=c("aveDnDs", "Expected"),
                        ylabel="dN/dS", title="dN/dS")

#+ fig.width=20
plot_dnds_across_genome(cols=c("dNdSWeightBySubst", "dNdSWeightByReads", "Expected"),
                        ylabel="dN/dS", title="dN/dS")

#+ fig.width=20
plot_dnds_across_genome(cols=c("dNdSWeightByReadsNoLowSyn", "dNdSWeightByReadsNoLowSubst", "Expected"),
                        ylabel="dN/dS", title="dN/dS")

#+ fig.width=20
plot_dnds_across_genome(cols=c("dN_minus_dS", "ExpectedMinus"),
                        ylabel="dN-dS", title="dN-dS")

#+ fig.width=20
plot_dnds_across_genome(cols=c("multisiteAvedNdS", "multisiteAvedNdSWeightBySubst", "multisitedNdSWeightBySubst", "ExpectedMultisite"),
                        ylabel="Multisite dN/dS", title="Multisite dN/dS")

#+ fig.width=20
plot_dnds_across_genome(cols=c("multisitedNdSWeightBySubst", "ExpectedMultisite"),
                        ylabel="Multisite dN/dS", title="Multisite dN/dS")



#' **Plot the Actual Codon Coverage across genome**
#' 
#+ fig.width=20
ggplot(collate_dnds, aes(x=CodonSite, y=CodonDepth) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Sequences with Unambiguous CodonDepth") + 
  ggtitle("Codon Coverage Across Genome")

#' **Plot the nonsynonymous substitutions across phylogeny**
#' 
fullDat$ExpectedSubst <- fullDat$ExpectedNonSyn + fullDat$ExpectedSyn
fullDatSubst <- reshape2:::melt(data=fullDat[, c("CodonSite", "ExpectedNonSyn", "ExpectedSyn", "ExpectedSubst", "NonSyn", "Syn", "Subst")], 
                                id.vars="CodonSite",
                                measure.vars=c("ExpectedNonSyn", "ExpectedSyn", "ExpectedSubst", "NonSyn", "Syn", "Subst"),
                                variable.name="substSource", value.name="total")
summary(fullDatSubst)
head(fullDatSubst)

#+ fig.width=20
ggplot(fullDatSubst[fullDatSubst$substSource %in% c("NonSyn", "ExpectedNonSyn"),], aes(x=CodonSite, y=total, color=substSource) ) + 
  geom_line() + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("Nonsynonymous Substitutions") + 
  ggtitle("Nonsynonymous Substitutions Across Phylogeny")

#' **Plot the synonymous substitutions across phylogeny**
#' 
#+ fig.width=20
ggplot(fullDatSubst[fullDatSubst$substSource %in% c("Syn", "ExpectedSyn"),], aes(x=CodonSite, y=total, color=substSource) ) + 
  geom_line() + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("Synonymous Substitutions") + 
  ggtitle("Synonymous Substitutions Across Phylogeny")

#' **Plot the substitutions across phylogeny**
#' 
#+ fig.width=20
ggplot(fullDatSubst[fullDatSubst$substSource %in% c("Subst", "ExpectedSubst"),], aes(x=CodonSite, y=total, color=substSource) ) + 
  geom_line() + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("Substitutions") + 
  ggtitle("Substitutions Across Phylogeny")

#' **Plot the Windows across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=CodonSite, y=Windows) ) + 
  geom_line() + 
  geom_smooth() + 
  xlab("Codon Site Along Genome") + 
  ylab("Windows") + 
  ggtitle("Windows Across Genome")


#' **Plot the Expected Omega rate across the genome**
#' 
#+ fig.width=20
ggplot(expected_dnds, aes(x=CodonSite, y=Omega) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dn/dS Expected") + 
  ggtitle("Expected Selection Along Genome")



#'**Concordance**
#'

# Returns a table of Lin's concordance correlation values
print_table_corr <- function() {
  site_dnds_corr <- sapply(c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "dNdSWeightByReadsNoLowSyn", "dNdSWeightByReadsNoLowSubst"),
         function(col) {
           dnds_ccc <- epi.ccc(fullDat[, col], fullDat$Expected)
           return (dnds_ccc$rho.c$est)
         })
  
  multisite_dnds_corr <- sapply(c("multisiteAvedNdS", "multisiteAvedNdSWeightBySubst", "multisitedNdSWeightBySubst"),
                           function(col) {
                             dnds_ccc <- epi.ccc(fullDat[, col], fullDat$ExpectedMultisite)
                             return (dnds_ccc$rho.c$est)
                           })
  
  site_dn_minus_ds_corr <- sapply(c("dN_minus_dS", "dN_minus_dSWeightByReadsNoLowSubst", "dN_minus_dSWeightByReads"),
                                function(col) {
                                  dnds_ccc <- epi.ccc(fullDat[, col], fullDat$ExpectedMinus)
                                  return (dnds_ccc$rho.c$est)
                                })
  corr_vals <- data.frame(Concordance=c(site_dnds_corr, multisite_dnds_corr, site_dn_minus_ds_corr))
}

table_corr <- print_table_corr()
kable(table_corr, format="html", caption="Concordance Correlation")


# Plots Lin's concordance correlation values by Non-Window Sequence Divergence
# Plot count of dn/ds values by Sequence Divergence
plot_table_corr_by_div <- function(nowindow_leaf_diversity) {
  
  nowindow_leaf_diversity$DivergeBin <- cut(nowindow_leaf_diversity$Diverge, 
                                            breaks=seq(0, max(nowindow_leaf_diversity$Diverge, na.rm=TRUE)+0.01, 0.01),
                                            include.lowest = TRUE, right=FALSE)
  full_div_combo <- merge(x=fullDat, 
                          y=subset(nowindow_leaf_diversity, select=c(CodonSite, DivergeBin)), 
                          by="CodonSite", all=TRUE, suffixes=c(".WinAve", ".NoWin"))
  
  corr_vals <-ddply(.data=full_div_combo,
                    .variables="DivergeBin",
                     .fun=function(x) {
                        data.frame(                                               
                          aveDnDs=epi.ccc(x$aveDnDs, x$Expected)$rho.c$est,
                          dNdSWeightBySubst=epi.ccc(x$dNdSWeightBySubst, x$Expected)$rho.c$est,
                          dNdSWeightByReads=epi.ccc(x$dNdSWeightByReads, x$Expected)$rho.c$est,
                          dNdSWeightByReadsNoLowSyn=epi.ccc(x$dNdSWeightByReadsNoLowSyn, x$Expected)$rho.c$est,
                          dNdSWeightByReadsNoLowSubst=epi.ccc(x$dNdSWeightByReadsNoLowSubst, x$Expected)$rho.c$est,
                          multisiteAvedNdS=epi.ccc(x$multisiteAvedNdS, x$ExpectedMultisite)$rho.c$est,
                          multisiteAvedNdSWeightBySubst=epi.ccc(x$multisiteAvedNdSWeightBySubst, x$ExpectedMultisite)$rho.c$est,
                          multisitedNdSWeightBySubst=epi.ccc(x$multisitedNdSWeightBySubst, x$ExpectedMultisite)$rho.c$est,
                          dN_minus_dS=epi.ccc(x$dN_minus_dS, x$ExpectedMinus)$rho.c$est,
                          dN_minus_dSWeightByReadsNoLowSubst=epi.ccc(x$dN_minus_dSWeightByReadsNoLowSubst, x$ExpectedMinus)$rho.c$est,
                          dN_minus_dSWeightByReads=epi.ccc(x$dN_minus_dSWeightByReads, x$ExpectedMinus)$rho.c$est)
                  })
  
  corr_vals_melt <- reshape2:::melt(data=corr_vals, id.vars="DivergeBin", variable.name="dnds_calc", value.name="Concordance")
  corr_vals_melt$Concordance[is.na(corr_vals_melt$Concordance)] <-  0
  fig <- ggplot(corr_vals_melt, aes(x=DivergeBin, y=Concordance)) + 
    geom_point() + 
    geom_line() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_wrap(~dnds_calc)
  print(fig)
  
}

#' What is the Concordance By Full Population Diversity
#' =======================================================================================
#' 
#+ fig.width=15, fig.height=12
 plot_table_corr_by_div(full_popn_conserve_dat$CodonDat)



#'
#' What is the Concordance By original Read Diversity
#' =======================================================================================
#' 
#+ fig.width=15, fig.height=12
 plot_table_corr_by_div(orig_conserve_dat$CodonDat)

#'
#' What is the Concordance By Aligned Read Diversity
#' =======================================================================================
#' 
#+ fig.width=15, fig.height=12
 plot_table_corr_by_div(aln_conserve_dat$CodonDat)


#' What Happens to Sites with Expected dN/dS < 1?
#' ================================================================
#'  
sites_exp_purify <- expected_dnds[expected_dnds$Omega < 1, "CodonSite"]
summary(collate_dnds[collate_dnds$CodonSite %in% sites_exp_purify, ])

#' What Happens to Sites with Divergence < `r mean(full_popn_conserve_dat$CodonDat$Diverge)` with Expected dN/dS < 1?
#' ================================================================
#'  
sites_low_act_div_exp_purify <- intersect(full_popn_conserve_dat$CodonDat[full_popn_conserve_dat$CodonDat$Diverge < mean(full_popn_conserve_dat$CodonDat$Diverge), "CodonSite"],
                                             sites_exp_purify)
summary(collate_dnds[collate_dnds$CodonSite %in% sites_low_act_div_exp_purify, ])


#' What Happens to Sites with Expected dN/dS > 1?

#' ================================================================
#'  
sites_exp_diversify <- expected_dnds[expected_dnds$Omega > 1, "CodonSite"]
summary(collate_dnds[collate_dnds$CodonSite %in% sites_exp_diversify, ])

#' What Happens to Sites with Divergence < `r mean(full_popn_conserve_dat$CodonDat$Diverge)`  with Expected dN/dS > 1?
#' ================================================================
#'  
sites_low_act_div_exp_diversify <- intersect(full_popn_conserve_dat$CodonDat[full_popn_conserve_dat$CodonDat$Diverge < mean(full_popn_conserve_dat$CodonDat$Diverge), "CodonSite"],
                                                sites_exp_diversify)
summary(collate_dnds[collate_dnds$CodonSite %in% sites_low_act_div_exp_diversify, ])


