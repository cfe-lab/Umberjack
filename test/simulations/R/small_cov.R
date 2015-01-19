# Plot coverage of simulated test data
library(ggplot2)
library(knitr)
library(reshape2)
library(plyr)
library(scales)


PICARD_WGS_METRICS_SKIP <- 6
PICARD_WGS_METRICS_ROWS <- 1
PICARD_COV_HISTO_SKIP <- 10





# READ IN CONFIGS
R_CONFIG_FILENAME <- "./small_cov.config"
rconfig<-read.table(R_CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))
NUM_INDIV <- as.numeric(rconfig[rconfig$key=="NUM_INDIV",]$val)
NUM_CODON_SITES <- as.numeric(rconfig[rconfig$key=="NUM_CODON_SITES",]$val)
ART_FOLD_COVER <- as.numeric(rconfig[rconfig$key=="ART_FOLD_COVER",]$val)
orig_err_free_cov_tsv <- rconfig[rconfig$key=="ORIG_ERR_FREE_COV_TSV",]$val
aln_err_free_cov_tsv <- rconfig[rconfig$key=="ALN_ERR_FREE_COV_TSV",]$val
orig_err_free_wgs_metrics_file <- rconfig[rconfig$key=="ORIG_ERR_FREE_WGS_METRICS",]$val
aln_err_free_wgs_metrics_file <- rconfig[rconfig$key=="ALN_ERR_FREE_WGS_METRICS",]$val
orig_cov_tsv <- rconfig[rconfig$key=="ORIG_COV_TSV",]$val
aln_cov_tsv <- rconfig[rconfig$key=="ALN_COV_TSV",]$val
orig_wgs_metrics_file <- rconfig[rconfig$key=="ORIG_WGS_METRICS",]$val
aln_wgs_metrics_file <- rconfig[rconfig$key=="ALN_WGS_METRICS",]$val
print(rconfig)


#'  NUM_INDIV=`r NUM_INDIV`
#'  NUM_CODON_SITES=`r NUM_CODON_SITES`
#'  ART_FOLD_COVER=`r ART_FOLD_COVER`


POPN_BP <- NUM_INDIV * NUM_CODON_SITES * 3
GENOME_BP <- NUM_CODON_SITES * 3

plot_cov <- function(orig_cov_tsv, aln_cov_tsv, qual) {

  
  # samtools depth files do not report positions with zero coverage.  Add the missing ref/positions
  Ref <- data.frame(Ref=paste0("otu", seq(1, NUM_INDIV)))
  Pos <- data.frame(Pos=seq(1, GENOME_BP))
  full_table <- merge(Ref, Pos, by=NULL, all.x=TRUE, all.y=TRUE)
  summary(full_table)
  head(full_table)
  
  orig_cov <- read.table(orig_cov_tsv, header=FALSE)
  colnames(orig_cov) <- c("Ref", "Pos", "Cov")
  summary(orig_cov)
  orig_cov_full <- merge(x=orig_cov, y=full_table, all=TRUE)
  orig_cov_full$Cov[is.na(orig_cov_full$Cov)] <- 0
  summary(orig_cov_full)
  dim(orig_cov_full)
  
  # Overall sum of reads covering each Pos across all individuals
  orig_cov_overall <- aggregate(formula=Cov~Pos, data=orig_cov_full, FUN=sum)
  summary(orig_cov_overall)
  head(orig_cov_overall)
  
  # For each (Pos, Cov) combination, find number of individals with that combination
  orig_cov_numindiv <- count(df=orig_cov[orig_cov$Cov > 0,], vars=c("Pos", "Cov"))
  colnames(orig_cov_numindiv)[grep("freq", colnames(orig_cov_numindiv))] <- "NumIndiv"
  summary(orig_cov_numindiv)
  head(orig_cov_numindiv)
  
  aln_cov <- read.table(aln_cov_tsv, header=FALSE)
  colnames(aln_cov) <- c("Ref", "Pos", "Cov")
  summary(aln_cov)
  
  aln_cov_full <- merge(aln_cov, Pos, by="Pos", all=TRUE)
  summary(aln_cov_full)
  head(aln_cov_full)
  aln_cov_full$Cov[is.na(aln_cov_full$Cov)] <- 0
  aln_cov_full$Ref[is.na(aln_cov_full$Ref)] <- "consensus"
  summary(aln_cov_full)
  dim(aln_cov_full)
  
  fig <-ggplot(orig_cov_numindiv, aes(x=Pos, y=Cov)) + 
    geom_tile(aes(fill = NumIndiv),colour = "white") + 
    scale_fill_gradient(low = "steelblue", high = "red") + 
    ggtitle(paste0("Heatmap of Total Individuals With Coverage At Each Position, ", qual))
  print(fig)
  
  fig <- ggplot(orig_cov_overall, aes(x=Pos, y=Cov)) + 
    geom_line() +  
    ggtitle(paste0("Genome Coverage of Original Reads, ", qual))
  print(fig)
  
  fig <- ggplot(aln_cov, aes(x=Pos, y=Cov)) + 
    geom_line() +  
    ggtitle(paste0("Genome Coverage of Aligned Reads, ", qual))
  print(fig)  
  
  orig_metrics <- paste0("Genome Length (bp)=", GENOME_BP, 
                    " Population Total Bases (bp)=", POPN_BP, 
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
get_metrics <- function (orig_wgs_metrics_file) {
  wgs_metrics <- read.table(orig_wgs_metrics_file, skip=PICARD_WGS_METRICS_SKIP, header=TRUE, nrows=PICARD_WGS_METRICS_ROWS)
  summary(wgs_metrics)
  wgs_metrics <- as.vector(wgs_metrics)
  metrics <- paste0("\n% Bases with MapQ < 20=", 100*wgs_metrics["PCT_EXC_MAPQ"],
                    "\n% Bases with Q < 20=", 100*wgs_metrics["PCT_EXC_BASEQ"],
                    "\n% Bases Overlap =", 100*wgs_metrics["PCT_EXC_OVERLAP"])
  return (metrics)
}

#' Error Free Reads
#' ====================================
#' 

plot_cov(orig_err_free_cov_tsv, aln_err_free_cov_tsv, "Error Free")


#' **ORIGINAL ERROR-FREE READ METRICS:**  
#' **`r get_metrics(orig_err_free_wgs_metrics_file)`**
#' 

#' **ALIGNED ERROR-FREE READ METRICS:**  
#' **`r get_metrics(aln_err_free_wgs_metrics_file)`**
#' 

#' Typical Miseq Quality Reads
#' =======================================
#' 
#' Coverage for Bases with >20 Quality and Reads >20 Mapping quality.

plot_cov(orig_cov_tsv, aln_cov_tsv, "Typical Quality")
#' **ORIGINAL READ METRICS:**  
#' **`r get_metrics(orig_wgs_metrics_file)`**
#' 
#' **ALIGNED READ METRICS:**  
#' **`r get_metrics(aln_wgs_metrics_file)`**
#' 