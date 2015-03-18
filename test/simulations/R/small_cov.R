#+ setup, include=FALSE
library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, width=1200)



# Plot coverage of simulated test data
library(ggplot2)
library(knitr)
library(reshape2)
library(plyr)
library(scales)



PICARD_WGS_METRICS_SKIP <- 6
PICARD_WGS_METRICS_ROWS <- 1
PICARD_COV_HISTO_SKIP <- 10



#' READ IN CONFIGS
#' =============================
#' 
R_CONFIG_FILENAME <- "./small_cov.config"
rconfig<-read.table(R_CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

FULL_POPN_CONSERVE_CSV <- rconfig[rconfig$key=="FULL_POPN_CONSERVE_CSV",]$val
ART_FOLD_COVER <- as.numeric(rconfig[rconfig$key=="ART_FOLD_COVER",]$val)
ORIG_ERR_FREE_COV_TSV <- rconfig[rconfig$key=="ORIG_ERR_FREE_COV_TSV",]$val
ALN_ERR_FREE_COV_TSV <- rconfig[rconfig$key=="ALN_ERR_FREE_COV_TSV",]$val
ORIG_ERR_FREE_WGS_METRICS <- rconfig[rconfig$key=="ORIG_ERR_FREE_WGS_METRICS",]$val
ALN_ERR_FREE_WGS_METRICS <- rconfig[rconfig$key=="ALN_ERR_FREE_WGS_METRICS",]$val
ORIG_ERR_FREE_CONSERVE_CSV <- rconfig[rconfig$key=="ORIG_ERR_FREE_CONSERVE_CSV",]$val
ALN_ERR_FREE_CONSERVE_CSV <- rconfig[rconfig$key=="ALN_ERR_FREE_CONSERVE_CSV",]$val
ORIG_COV_TSV <- rconfig[rconfig$key=="ORIG_COV_TSV",]$val
ALN_COV_TSV <- rconfig[rconfig$key=="ALN_COV_TSV",]$val
ORIG_WGS_METRICS <- rconfig[rconfig$key=="ORIG_WGS_METRICS",]$val
ALN_WGS_METRICS <- rconfig[rconfig$key=="ALN_WGS_METRICS",]$val
ORIG_CONSERVE_CSV <- rconfig[rconfig$key=="ORIG_CONSERVE_CSV",]$val
ALN_CONSERVE_CSV <- rconfig[rconfig$key=="ALN_CONSERVE_CSV",]$val
INDELIBLE_RATES_CSV <- rconfig[rconfig$key=="INDELIBLE_RATES_CSV", ]$val
CMP_READ_ERR_CSV <- rconfig[rconfig$key=="CMP_READ_ERR_CSV", ]$val

#+ results='asis'
kable(rconfig, format="html", caption="Configs")


#' Read in Diversity CSVs
#' ================================
#' 

#+
# Takes a *.conserve.csv file with nucleotide stats for unsliced sequences
# Expects cols: "NucSite", "Conserve", "Entropy", "NucDepth", "CodonDepth"  
# Returns a dataframe with elements:
# - Conserve, Diverge, Entropy, CodonDepth, NucDepth
fill_nuc_conserve <- function(nuc_conserve_csv) {
  # Per-nucleotide conservation, Divergence, Entropy
  conserve_nuc <- read.table(nuc_conserve_csv, header=TRUE, sep=",")
  conserve_nuc$Diverge <- 1 - conserve_nuc$Conserve
  return (conserve_nuc)
}

orig_conserve <- fill_nuc_conserve(ORIG_CONSERVE_CSV)
summary(orig_conserve)

orig_errfree_conserve <- fill_nuc_conserve(ORIG_ERR_FREE_CONSERVE_CSV)
summary(orig_errfree_conserve)

aln_conserve <- fill_nuc_conserve(ALN_CONSERVE_CSV)
summary(aln_conserve)

aln_errfree_conserve <- fill_nuc_conserve(ALN_ERR_FREE_CONSERVE_CSV)
summary(aln_errfree_conserve)

full_popn_conserve <- fill_nuc_conserve(FULL_POPN_CONSERVE_CSV)
summary(full_popn_conserve)

# Find the total individuals and total codon sites from the full population conservation csv
NUM_NUC_SITES <- max(full_popn_conserve$NucSite)
NUM_CODON_SITES <- ceiling(NUM_NUC_SITES / 3)
NUM_INDIV <- max(full_popn_conserve$NucDepth)
POPN_BP <- NUM_INDIV * NUM_NUC_SITES

#'  NUM_NUC_SITES=`r NUM_NUC_SITES`
#'  NUM_CODON_SITES=`r NUM_CODON_SITES`
#'  NUM_INDIV=`r NUM_INDIV`
#'  POPN_BP=`r POPN_BP`
#'  

indelible <- read.table(INDELIBLE_RATES_CSV, header=TRUE, sep=",")
summary(indelible)

#' Read in Sequencing Error CSVs
#' ================================
#' 


# Get merged read stats on sequencing error, length, N's
read_stats <- read.table(CMP_READ_ERR_CSV, header=TRUE, sep=",")
read_stats$Len <- read_stats$End - read_stats$Start + 1
summary(read_stats)
summary(read_stats[read_stats$Source == "Orig",])
summary(read_stats[read_stats$Source == "OrigErrFree",])
summary(read_stats[read_stats$Source == "Aln",])
summary(read_stats[read_stats$Source == "AlnErrFree",])




#' Function Definitions
#' ==============================
#' 
 
# Boxplot the actual divergence, entropy at each Indelible mutation scaling rate
plot_act_diversity_by_indelible_mut <- function(full_popn_conserve, orig_conserve, orig_errfree_conserve,
                                                aln_conserve, aln_errfree_conserve, indelible) {
  all <- rbind(data.frame(Source="Full", subset(full_popn_conserve, select=c("NucSite", "Diverge", "Entropy"))),
               data.frame(Source="Orig", subset(orig_conserve, select=c("NucSite", "Diverge", "Entropy"))),
               data.frame(Source="OrigErrFree", subset(orig_errfree_conserve, select=c("NucSite", "Diverge", "Entropy"))),
               data.frame(Source="Aln", subset(aln_conserve, select=c("NucSite", "Diverge", "Entropy"))),
               data.frame(Source="AlnErrFree", subset(aln_errfree_conserve, select=c("NucSite", "Diverge", "Entropy"))))
  
  all$CodonSite <- ceiling(all$NucSite / 3)
  all <- merge(x=all, y=indelible, by.x="CodonSite", by.y="Site", all.x=TRUE, all.y=FALSE)
  
  
  fig <- ggplot(all, aes(x=as.factor(Scaling_factor), y=Diverge, color=Source)) +
    geom_boxplot() + 
    xlab("\n Indelible Mutation Scaling Factor") + 
    ylab("Nucleotide Site Divergence from Consensus\n") + 
    ggtitle("Boxplot Actual Divergence by Indelible Mutation Scaling Factor")
  print(fig)
  
  
  fig <- ggplot(all, aes(x=as.factor(Scaling_factor), y=log10(Entropy), color=Source)) +
    geom_boxplot() + 
    xlab("\n Indelible Mutation Scaling Factor") + 
    ylab("Log10 Nucleotide Site Metric Entropy\n") + 
    ggtitle("Boxplot Actual Metric Entropy by Indelible Mutation Scaling Factor")
  print(fig)
}

# Boxplot the coverage at each Indelible mutation scaling rate
plot_coverage_by_indelible_mut <- function(orig_read_cov_tsv=ORIG_COV_TSV, 
                                           aln_read_cov_tsv=ALN_COV_TSV, 
                                           orig_read_errfree_cov_tsv=ORIG_ERR_FREE_COV_TSV, 
                                           aln_read_err_free_cov_tsv=ALN_ERR_FREE_COV_TSV,
                                           indelible) {
  
  orig_cov <- read.table(orig_read_cov_tsv, header=FALSE, col.names=c("Ref", "NucSite", "Cov"))
  summary(orig_cov)
  
  # Overall sum of reads covering each Pos across all individuals
  orig_cov_full <- aggregate(formula=Cov~NucSite, data=orig_cov, FUN=sum)
  orig_cov_full <- merge(orig_cov_full, data.frame(NucSite=c(1:NUM_NUC_SITES)), by="NucSite", all=TRUE)
  orig_cov_full$Cov[is.na(orig_cov_full$Cov)] <- 0
  print(summary(orig_cov_full))
  dim(orig_cov_full)
  
  orig_cov_errfree <- read.table(orig_read_errfree_cov_tsv, header=FALSE, col.names=c("Ref", "NucSite", "Cov"))
  orig_cov_errfree_full <- aggregate(formula=Cov~NucSite, data=orig_cov_errfree, FUN=sum)
  orig_cov_errfree_full <- merge(orig_cov_full, data.frame(NucSite=c(1:NUM_NUC_SITES)), by="NucSite", all=TRUE)
  orig_cov_errfree_full$Cov[is.na(orig_cov_errfree_full$Cov)] <- 0
  print(summary(orig_cov_errfree_full))
  dim(orig_cov_errfree_full)
  
  aln_cov <- read.table(aln_read_cov_tsv, header=FALSE, col.names=c("Ref", "NucSite", "Cov"))
  summary(aln_cov)
  # samtools depth files do not report positions with zero coverage.  Add the missing ref/positions
  aln_cov_full <- merge(aln_cov, data.frame(NucSite=c(1:NUM_NUC_SITES)), by="NucSite", all=TRUE)
  summary(aln_cov_full)
  head(aln_cov_full)
  aln_cov_full$Cov[is.na(aln_cov_full$Cov)] <- 0
  print(summary(aln_cov_full))
  dim(aln_cov_full)
  
  aln_cov_errfree <- read.table(aln_read_err_free_cov_tsv, header=FALSE, col.names=c("Ref", "NucSite", "Cov"))
  summary(aln_cov_errfree)
  # samtools depth files do not report positions with zero coverage.  Add the missing ref/positions
  aln_cov_errfree_full <- merge(aln_cov_errfree, data.frame(NucSite=c(1:NUM_NUC_SITES)), by="NucSite", all=TRUE)
  summary(aln_cov_errfree_full)
  head(aln_cov_errfree_full)
  aln_cov_errfree_full$Cov[is.na(aln_cov_errfree_full$Cov)] <- 0
  print(summary(aln_cov_errfree_full))
  dim(aln_cov_errfree_full)
  
  all <- rbind(data.frame(Source="Orig", subset(orig_cov_full, select=c("NucSite", "Cov"))),
               data.frame(Source="Aln", subset(orig_cov_errfree_full, select=c("NucSite", "Cov"))),
               data.frame(Source="OrigErrFree", subset(aln_cov_full, select=c("NucSite", "Cov"))),
               data.frame(Source="AlnErrFree", subset(aln_cov_errfree_full, select=c("NucSite", "Cov"))))
  all$CodonSite <- ceiling(all$NucSite / 3)
  all <- merge(x=all, y=indelible, by.x="CodonSite", by.y="Site", all.x=TRUE, all.y=FALSE)
  print(summary(all))
  print(dim(all))
  fig <- ggplot(all, aes(x=as.factor(Scaling_factor), y=Cov, color=Source)) +
    geom_boxplot() + 
    xlab("\n Indelible Mutation Scaling Factor") + 
    ylab("Read Coverage\n") + 
    ggtitle("Boxplot Read Coverage Indelible Mutation Scaling Factor")
  print(fig)
}

# Plots divergence, entropy at each nucleotide position, codon position
plot_conserve <- function(conserve_csv, title) {
  # Col: "NucSite", "Conserve", "Entropy", "NucDepth", "CodonDepth"
  conserve <- read.table(conserve_csv, header=TRUE, sep=",")
  conserve$Diverge <- 1 - conserve$Conserve
  print(summary(conserve))
  
  fig <- ggplot(conserve, aes(x=NucSite, y=Diverge)) + 
    geom_line() + 
    geom_smooth() + 
    xlab("\n Nucleotide Site") + 
    ylab("Fraction Divergence\n") + 
    ggtitle(paste0(title, " Divergence At Each Nucleotide Site"))
  print(fig)
  
  fig <- ggplot(conserve, aes(x=NucSite, y=Entropy)) + 
    geom_line() + 
    geom_smooth() + 
    xlab("\n Nucleotide Site") + 
    ylab("Metric Entropy\n") + 
    ggtitle(paste0(title, " Metric Entropy At Each Nucleotide Site"))
  print(fig)
  
  conserve_codon <- data.frame(CodonSite=c(1:NUM_CODON_SITES))  # 1-based  codon sites
  
  conserve_codon$Diverge <- sapply(conserve_codon$CodonSite, function(codonsite) {
    start_nuc_site <- (codonsite * 3) -2  # 1-based nucleotide sites
    end_nuc_site <- codonsite * 3
    mean(conserve[conserve$NucSite >= start_nuc_site & conserve$NucSite <= end_nuc_site, "Diverge"])
  })
  
  conserve_codon$Entropy <- sapply(conserve_codon$CodonSite, function(codonsite) {
    start_nuc_site <- (codonsite * 3) -2  # 1-based nucleotide sites
    end_nuc_site <- codonsite * 3
    mean(conserve[conserve$NucSite >= start_nuc_site & conserve$NucSite <= end_nuc_site, "Entropy"])
  })
  print(summary(conserve_codon))
  
  fig <- ggplot(conserve_codon, aes(x=CodonSite, y=Diverge)) + 
    geom_line() + 
    geom_smooth() + 
    xlab("\n Codon Site") + 
    ylab("Fraction Divergence\n") + 
    ggtitle(paste0(title, " Divergence At Each Codon Site"))
  print(fig)
  
  fig <- ggplot(conserve_codon, aes(x=CodonSite, y=Entropy)) + 
    geom_line() + 
    geom_smooth() + 
    xlab("\n Codon Site") + 
    ylab("Metric Entropy\n") + 
    ggtitle(paste0(title, " Metric Entropy At Each Codon Site"))
  print(fig)
}

# Compares on the same graph, the diversity, entropy at nucleotide level and codon level between alignments for 
#   full population, original reads, aligned reads
plot_cmp_conserve <- function(full_popn_conserve_csv, orig_read_conserve_csv, aln_read_conserve_csv, 
                              orig_read_errfree_conserve_csv, aln_read_err_free_conserve_csv) {  
  
  full_popn_conserve <- read.table(full_popn_conserve_csv, header=TRUE, sep=",")
  full_popn_conserve$Diverge <- 1- full_popn_conserve$Conserve
  
  orig_read_conserve <- read.table(orig_read_conserve_csv, header=TRUE, sep=",")
  orig_read_conserve$Diverge <- 1- orig_read_conserve$Conserve
  
  orig_read_err_free_conserve <- read.table(orig_read_errfree_conserve_csv, header=TRUE, sep=",")
  orig_read_err_free_conserve$Diverge <- 1- orig_read_err_free_conserve$Conserve
  
  aln_read_conserve <- read.table(aln_read_conserve_csv, header=TRUE, sep=",")
  aln_read_conserve$Diverge <- 1- aln_read_conserve$Conserve
  
  aln_read_err_free_conserve <- read.table(aln_read_err_free_conserve_csv, header=TRUE, sep=",")
  aln_read_err_free_conserve$Diverge <- 1- aln_read_err_free_conserve$Conserve
  
  # Plot full popn, orig reads, aln reads on same graph
  # Cols:  Source, NucSite, Diverge, Conserve, Entropy
  all_conserve <- rbind(data.frame(Source="Full", full_popn_conserve),
                        data.frame(Source="Orig", orig_read_conserve),
                        data.frame(Source="Aln", aln_read_conserve),
                        data.frame(Source="OrigErrFree", orig_read_err_free_conserve),
                        data.frame(Source="AlnErrFree", aln_read_err_free_conserve))
  summary(all_conserve)
  
  
  # Plot Diversity Along Genome
  sapply(c("Diverge", "Entropy", "NucDepth", "CodonDepth"), function(colname) {
    fig <- ggplot(all_conserve, aes(x=NucSite, color=Source, lty=Source)) + 
      geom_smooth(aes_string(y=colname)) + 
      xlab("\n Nucleotide Site") + 
      ylab(paste0(colname, "\n")) + 
      ggtitle(paste0("Compare Smoothed ", colname, " At Each Nucleotide Site"))
    print(fig)
  })
  
  
  # Plot Density of Sites For Diversity
  sapply(c("Diverge", "Entropy", "NucDepth", "CodonDepth"), function(colname) {
    fig <- ggplot(all_conserve, aes(x=NucSite, color=Source)) + 
      geom_density(aes_string(x=colname)) + 
      xlab(paste0("\n", colname)) + 
      ylab(paste0("Density of Sites\n")) + 
      ggtitle(paste0("Compare Density of Sites For ", colname))
    print(fig)
  })
  
  # Plot Cumulative Frequency of Sites For Diversity
  sapply(c("Diverge", "Entropy", "NucDepth", "CodonDepth"), function(colname) {
    fig <- ggplot(all_conserve, aes_string(x=colname, color="Source", lty="Source")) + 
      stat_ecdf() + 
      xlab(paste0("\n", colname)) + 
      ylab(paste0("Cumulative Distribution of Sites\n")) + 
      ggtitle(paste0("Compare Cumulative Distribution of Sites For ", colname))
    print(fig)
  })
  
  # Scatter Orig Read Diversity Vs Full Population Diversity    
  full_popn_melt <- reshape2:::melt(data=subset(full_popn_conserve, select=c(NucSite, CodonDepth, Diverge, Entropy)),
                                         id.vars="NucSite",
                                         measure.vars=c("Diverge", "Entropy"),
                                         variable.name="DiversityVar.Full", 
                                         value.name="DiversityVal.Full")
  full_popn_melt$DiversityVar.Full <- paste0(full_popn_melt$DiversityVar.Full, ".Full")
  
  orig_melt <-  rbind(data.frame(Source="Orig", subset(orig_read_conserve, select=c(NucSite, CodonDepth, Diverge, Entropy))),
                      data.frame(Source="OrigErrFree", subset(orig_read_err_free_conserve, select=c(NucSite, CodonDepth, Diverge, Entropy))))
  orig_melt <- reshape2:::melt(data=orig_melt,
                                    id.vars=c("Source", "NucSite"),
                                    measure.vars=c("CodonDepth", "Diverge", "Entropy"),
                                    variable.name="DiversityVar.Orig", 
                                    value.name="DiversityVal.Orig")
  orig_melt$DiversityVar.Orig <- paste0(orig_melt$DiversityVar.Orig, ".Orig")
  
  full_read_combo <- merge(x=full_popn_melt, y=orig_melt, by="NucSite", all=TRUE)
  
  fig <- ggplot(full_read_combo, aes(x=DiversityVal.Full, y=DiversityVal.Orig, color=Source)) + 
    facet_grid(DiversityVar.Orig~DiversityVar.Full, scales="free") + 
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm", aes(linetype=Source)) + 
    xlab("\nFull Population Sequence Diversity") + 
    ylab("Original Reads Sequence Diversity\n") + 
    ggtitle("Scatterplot Original Reads Sequence Diversity Vs Full Population Sequence Diversity")
  print(fig)
  
  aln_melt <-  rbind(data.frame(Source="Aln", subset(aln_read_conserve, select=c(NucSite, CodonDepth, Diverge, Entropy))),
                      data.frame(Source="AlnErrFree", subset(aln_read_err_free_conserve, select=c(NucSite, CodonDepth, Diverge, Entropy))))
  aln_melt <- reshape2:::melt(data=aln_melt,
                               id.vars=c("Source", "NucSite"),
                               measure.vars=c("CodonDepth", "Diverge", "Entropy"),
                               variable.name="DiversityVar.Aln", 
                               value.name="DiversityVal.Aln")
  aln_melt$DiversityVar.Aln <- paste0(aln_melt$DiversityVar.Aln, ".Aln")
  full_read_combo <- merge(x=full_popn_melt, y=aln_melt, by="NucSite", all=TRUE)
  
  fig <- ggplot(full_read_combo, aes(x=DiversityVal.Full, y=DiversityVal.Aln, color=Source)) + 
    facet_grid(DiversityVar.Aln~DiversityVar.Full, scales="free") + 
    geom_point(shape=1, alpha=0.1, aes(linetype=Source)) + 
    geom_smooth(method="lm") + 
    xlab("\nFull Population Sequence Diversity") + 
    ylab("Aligned Reads Sequence Diversity\n") + 
    ggtitle("Scatterplot Aligned Reads Sequence Diversity Vs Full Population Sequence Diversity")
  print(fig)
  
}
  




plot_cov <- function(ORIG_COV_TSV, ALN_COV_TSV, qual) {
  
  orig_cov <- read.table(ORIG_COV_TSV, header=FALSE, col.names=c("Ref", "NucSite", "Cov"))
  summary(orig_cov)
  
  # samtools depth files do not report positions with zero coverage.  Add the missing ref/positions
  # We need to get stats that include the positions with zero coverage.
  Ref <- data.frame(Ref=paste0("otu", seq(1, NUM_INDIV)))
  NucSite <- data.frame(NucSite=seq(1, NUM_NUC_SITES))
  full_table <- merge(Ref, NucSite, by=NULL, all.x=TRUE, all.y=TRUE)
  summary(full_table)
  head(full_table)
  orig_cov_full <- merge(x=orig_cov, y=full_table, all=TRUE)
  orig_cov_full$Cov[is.na(orig_cov_full$Cov)] <- 0
  print(summary(orig_cov_full))
  dim(orig_cov_full)

  
  # Overall sum of reads covering each Pos across all individuals
  orig_cov_overall <- aggregate(formula=Cov~NucSite, data=orig_cov, FUN=sum)
  summary(orig_cov_overall)
  head(orig_cov_overall)
  
  aln_cov <- read.table(ALN_COV_TSV, header=FALSE, col.names=c("Ref", "NucSite", "Cov"))
  summary(aln_cov)
  
  # samtools depth files do not report positions with zero coverage.  Add the missing ref/positions
  aln_cov_full <- merge(aln_cov, data.frame(NucSite=c(1:NUM_NUC_SITES)), by="NucSite", all=TRUE)
  summary(aln_cov_full)
  head(aln_cov_full)
  aln_cov_full$Cov[is.na(aln_cov_full$Cov)] <- 0
  print(summary(aln_cov_full))
  dim(aln_cov_full)

  
  orig_metrics <- paste0(
                    " Original Average Coverage Per Individual= ", mean(orig_cov_full$Cov),
                    " Original Stddev Coverage Per Individual= ", sd(orig_cov_full$Cov),
                    " Original Median Coverage Per Individual= ", median(orig_cov_full$Cov)
                    )
  print(orig_metrics)
  
  aln_metrics <- paste0( "Genome Aligned Average Coverage= ", mean(aln_cov_full$Cov),
                         " Genome Aligned Stddev Coverage= ", sd(aln_cov_full$Cov),
                         " Genome Aligned Median Coverage= ", median(aln_cov_full$Cov)
  )
  print (aln_metrics)
}

# NB:  For Picard WGS metrics, overlaps are not included in coverage, neither are bases with Q<20, bases with mapq < 20
get_metrics <- function (ORIG_WGS_METRICS) {
  wgs_metrics <- read.table(ORIG_WGS_METRICS, skip=PICARD_WGS_METRICS_SKIP, header=TRUE, nrows=PICARD_WGS_METRICS_ROWS)
  summary(wgs_metrics)
  wgs_metrics <- as.vector(wgs_metrics)
  metrics <- paste0("\n% Bases with MapQ < 20=", 100*wgs_metrics["PCT_EXC_MAPQ"],
                    "\n% Bases with Q < 20=", 100*wgs_metrics["PCT_EXC_BASEQ"],
                    "\n% Bases Overlap =", 100*wgs_metrics["PCT_EXC_OVERLAP"])
  return (metrics)
}

# This takes a loooong time due to the aggregation of the massive original read coverage dataframe.
# Compare Coverage of Full population, original reads (typical, err free), aligned reads (typical, error free) on same graph
plot_cmp_cov <- function(orig_read_cov_tsv, aln_read_cov_tsv, orig_read_errfree_cov_tsv, aln_read_err_free_cov_tsv) {
  
  orig_read_cov <- read.table(orig_read_cov_tsv, header=FALSE, sep="\t", col.names=c("Ref", "NucSite", "Cov"))
  orig_read_cov <- aggregate(formula=Cov~NucSite, data=orig_read_cov, FUN=sum)
  
  orig_read_errfree_cov <- read.table(orig_read_errfree_cov_tsv, header=FALSE, sep="\t", col.names=c("Ref", "NucSite", "Cov"))
  orig_read_errfree_cov <- aggregate(formula=Cov~NucSite, data=orig_read_errfree_cov, FUN=sum)
  
  aln_read_cov <- read.table(aln_read_cov_tsv, header=FALSE, sep="\t", col.names=c("Ref", "NucSite", "Cov"))
  aln_read_errfree_cov <- read.table(aln_read_err_free_cov_tsv, header=FALSE, sep="\t", col.names=c("Ref", "NucSite", "Cov"))
  
  # Plot per-population coverage along genome
  all_cov <- rbind(data.frame(Source="Orig", subset(orig_read_cov, select=c(NucSite, Cov))),
                   data.frame(Source="Aln", subset(aln_read_cov, select=c(NucSite, Cov))),
                   data.frame(Source="OrigErrFree", subset(orig_read_errfree_cov, select=c(NucSite, Cov))),
                   data.frame(Source="AlnErrFree", subset(aln_read_errfree_cov, select=c(NucSite, Cov)))
                   )
  summary(all_cov)

  fig <- ggplot(all_cov, aes(x=NucSite, y=Cov, color=Source, linetype=Source)) + 
    geom_line() +  
    xlab("\n Nucleotide Site") + 
    xlab("\n Per-Population Fold Coverage") + 
    ggtitle("Compare Per-Population Coverage Along Genome")
  print(fig)

  
}

#'
#' Do Nucleotide Sites Obey the Indelible Mutation Scaling Factors?
#' =====================================================================
#'  
#'  
plot_act_diversity_by_indelible_mut(full_popn_conserve, orig_conserve, orig_errfree_conserve,
                                    aln_conserve, aln_errfree_conserve, indelible)


#'
#' Does Indelible Mutation Scaling Factors Affect Read Coverage?
#' =====================================================================
#'  
#'  
plot_coverage_by_indelible_mut(orig_read_cov_tsv=ORIG_COV_TSV, 
                               aln_read_cov_tsv=ALN_COV_TSV, 
                               orig_read_errfree_cov_tsv=ORIG_ERR_FREE_COV_TSV, 
                               aln_read_err_free_cov_tsv=ALN_ERR_FREE_COV_TSV,
                               indelible)


#' Compare Full Population, Original Reads (Typical, Error Free), Aligned Reads (Typical, Error Free)
#' ====================================
#' 
#+ fig.width=20, warning=FALSE
plot_cmp_conserve(full_popn_conserve_csv=FULL_POPN_CONSERVE_CSV, 
                  orig_read_conserve_csv=ORIG_CONSERVE_CSV, 
                  aln_read_conserve_csv=ALN_CONSERVE_CSV, 
                  orig_read_errfree_conserve_csv=ORIG_ERR_FREE_CONSERVE_CSV, 
                  aln_read_err_free_conserve_csv=ALN_ERR_FREE_CONSERVE_CSV)

  
#+ fig.width=20, warning=FALSE
plot_cmp_cov(orig_read_cov_tsv=ORIG_COV_TSV, 
             aln_read_cov_tsv=ALN_COV_TSV, 
             orig_read_errfree_cov_tsv=ORIG_ERR_FREE_COV_TSV, 
             aln_read_err_free_cov_tsv=ALN_ERR_FREE_COV_TSV)

#' Full Population
#' ====================================
#' 
#'  Divergence, Entropy Across Genome for Full Population (No Sequencing)
#'  
#'  
plot_conserve(FULL_POPN_CONSERVE_CSV, "Full Population")




#' Error Free Reads
#' ====================================
#' 
#+ fig.width=20
plot_cov(ORIG_ERR_FREE_COV_TSV, ALN_ERR_FREE_COV_TSV, "Error Free")

#' **ORIGINAL ERROR-FREE READ METRICS:**  
#' **`r get_metrics(ORIG_ERR_FREE_WGS_METRICS)`**
#' 

#' **ALIGNED ERROR-FREE READ METRICS:**  
#' **`r get_metrics(ALN_ERR_FREE_WGS_METRICS)`**
#' 

plot_conserve(ORIG_ERR_FREE_CONSERVE_CSV, "Original Error Free Reads")
plot_conserve(ALN_ERR_FREE_CONSERVE_CSV, "Aligned Error Free Reads")

#' Typical Miseq Quality Reads
#' =======================================
#' 
#' Coverage for Bases with >20 Quality and Reads >20 Mapping quality.

plot_cov(ORIG_COV_TSV, ALN_COV_TSV, "Typical Quality")
#' **ORIGINAL READ METRICS:**  
#' **`r get_metrics(ORIG_WGS_METRICS)`**
#' 
#' **ALIGNED READ METRICS:**  
#' **`r get_metrics(ALN_WGS_METRICS)`**
#' 

plot_conserve(ORIG_CONSERVE_CSV, "Original Reads")
plot_conserve(ALN_CONSERVE_CSV, "Aligned Reads")
