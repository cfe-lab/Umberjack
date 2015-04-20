library(knitr)
opts_chunk$set(progress=FALSE, verbose=FALSE, warning=FALSE, message=FALSE, width=1500)

library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)
library(plyr)



# Read locations of input files from local umberjack_unit_test.config file
CONFIG_FILENAME <- "./umberjack_unit_test.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

actual_dnds_filename <- config[config$key=="ACTUAL_DNDS_FILENAME",]$val
expected_dnds_filename <-  config[config$key=="EXPECTED_DNDS_FILENAME",]$val



indelible_dnds_filename <- config[config$key=="INDELIBLE_DNDS_FILENAME", ]$val

#'  ACTUAL_DNDS_FILENAME=`r actual_dnds_filename`
#'  
#'  EXPECTED_DNDS_FILENAME=`r expected_dnds_filename`
#'  
#'  INDELIBLE_DNDS_FILENAME=`r indelible_dnds_filename`
#'  
actual_dnds_file <- file(actual_dnds_filename, open="rt")
comments <- readLines(actual_dnds_file, 1) # Read one line 
close(actual_dnds_file)

args <- unlist(strsplit(comments, ','))




#' *`r comments`**
#' -----------------------------
#' 
#' 
actual_dnds <- read.table(actual_dnds_filename, header=TRUE, na.strings="None", comment.char = "#", sep=",")
dim(actual_dnds)
head(actual_dnds)
str(actual_dnds)
summary(actual_dnds)

# Cols: Observed S Changes  Observed NS Changes	E[S Sites]	E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS	dn/ds
# Parse the expected dnds filename for it start and end nucleotide positions (1-based)

expected_dnds <- read.table(expected_dnds_filename, header=TRUE, sep="\t")  # Site  Interval	Scaling_factor	Rate_class	Omega
expected_dnds$Site <- as.numeric(rownames(expected_dnds))
expected_dnds$Omega <- expected_dnds$dN/expected_dnds$dS
expected_dnds$Omega[expected_dnds$dS == 0] <- NA
expected_dnds$Subs <- expected_dnds$Observed.S.Changes + expected_dnds$Observed.NS.Changes
dim(expected_dnds)
head(expected_dnds)
str(expected_dnds)
summary(expected_dnds)


# check consistency
if (!all(expected_dnds$Site == actual_dnds$Site)) {
  stop(paste0("Actual dnds csv ", actual_dnds_filename, " and Expected dnds csv ", expected_dnds_filename, " sites not in sync"))
}


#' Compare Full Population INDELible dN/dS with Full Population HyPhy dN/dS
#' ============================================================================

# Read in Indelible's intended dN/dS
# Cols: Site,Interval,Scaling_factor,Rate_class,Omega.  These columns are not ordered by site.
indelible_dnds <- read.table(indelible_dnds_filename, header=TRUE, sep=",")
indelible_dnds$Scaling_factor <- as.factor(indelible_dnds$Scaling_factor)
summary(indelible_dnds)
head(indelible_dnds)

#' **Plot indelible's intended dN/dS vs the actual dN/dS as determined by HyPhy**
#' 
indelhyphy <- merge(x=indelible_dnds, y=expected_dnds[, c("Site", "Subs", "Omega")], 
                    by.x="Site", by.y="Site", suffixes=c(".indelible", ".hyphy"), all=FALSE)
colnames(indelhyphy)[grep("Subs", colnames(indelhyphy))] <- "Subs.hyphy"
summary(indelhyphy)
indelhyphy <- na.omit(indelhyphy)

#' **Scatterplot Full Population HyPhy Actual dN/dS vs Indelible Expected dN/dS**
ggplot(indelhyphy, aes(x=Omega.indelible, y=Omega.hyphy)) + 
  geom_point(shape=1, alpha=0.5, na.rm=TRUE) +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
  geom_abline(slope=1) + 
  ylab("HyPhy dN/dS\n") + 
  xlab("\nIndelible dN/dS") + 
  #coord_fixed(ratio=1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  ggtitle("HyPhy dN/dS Vs Indelible dN/dS on Full Population")

#' **Boxplot Full Population HyPhy Actual dN/dS vs Indelible Expected dN/dS**
ggplot(indelhyphy, aes(x=as.factor(Omega.indelible), y=Omega.hyphy)) + 
  geom_boxplot() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
  ylab("HyPhy dN/dS\n") + 
  xlab("\nIndelible dN/dS") + 
  scale_y_continuous(breaks=seq(0, max(indelhyphy$Omega.hyphy), 0.5)) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=14) ,
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("HyPhy dN/dS Vs Indelible dN/dS on Full Population")


#' **Boxplot Full Population HyPhy Actual Subs vs Indelible Mutation Scaling Rate**
ggplot(indelhyphy, aes(x=Scaling_factor, y=Subs.hyphy)) + 
  geom_boxplot() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
  ylab("HyPhy Substitutions\n") + 
  xlab("\nIndelible Mutation Scaling Rate") + 
  ggtitle("HyPhy Subs Vs Indelible Subs on Full Population")



#' Compare Codon Site dN/dS Averaged Over Windows With Full Population HyPhy dN/dS
#' ============================================================================
#' 
#' **Scatterplot actual vs expected dn ds together**

fullDat <- merge(x=actual_dnds, 
                 y=expected_dnds[, c("Site", "Omega", "Scaled.dN.dS", "Observed.S.Changes", "Observed.NS.Changes")], 
                 by="Site", all=FALSE, sort=TRUE)
fullDat <- merge(x=fullDat, 
                 y=indelible_dnds[, c("Site", "Scaling_factor")], 
                 by="Site", all.x=TRUE, all.y=FALSE, sort=TRUE)
colnames(fullDat)[grep("Omega", colnames(fullDat))] <- "Expected"
colnames(fullDat)[grep("Scaled.dN.dS", colnames(fullDat))] <- "ExpectedMinus"
colnames(fullDat)[grep("Observed.S.Changes", colnames(fullDat))] <- "ExpectedSyn"
colnames(fullDat)[grep("Observed.NS.Changes", colnames(fullDat))] <- "ExpectedNonSyn"
summary(fullDat)

#' **Smoothed Scatterplot of Site dn/ds across the genome**
#' 
scatterplot_actual_v_expected <- function(expected_colname, expected_title, actual_colname, actual_title) {
  fig <- ggplot(fullDat  , aes_string(x=expected_colname, y=actual_colname)) + 
    geom_point(alpha=0.5, na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Inferred ", actual_title, "\n Vs Expected ", expected_title) )
  print(fig)
  
  actual_col <- fullDat[, actual_colname]
  ave <- mean(actual_col, na.rm=TRUE) 
  most_outlier <- max(abs( c(max(actual_col - ave, na.rm=TRUE), 
                             min(actual_col - ave, na.rm=TRUE))))
  stdev <- sd(actual_col, na.rm=TRUE)
  if (most_outlier > 4 * stdev) {
    fig <- ggplot(fullDat[abs(actual_col - ave) < 3 * stdev, ]  , aes_string(x=expected_colname, y=actual_colname)) + 
      geom_point(alpha=0.5, na.rm=TRUE) +
      geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
      geom_abline(slope=1) + 
      ylab("Inferred\n") + 
      xlab("\nExpected") + 
      #coord_fixed(ratio=1) + 
      theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
      ggtitle(paste0("No Outliers: Inferred ", actual_title, "\n Vs Expected ", expected_title) )
    print(fig)
  }
}
#Cols: Ref  Site  Windows  Codons  NonSyn  Syn	Subs	dNdSWeightByReads	dNdSWeightBySubs	dNdSWeightByReadsNoLowSub	dNdSWeightBySubsNoLowSub	dnMinusDsWeightByReads	dnMinusDsWeightBySubs	dnMinusDsWeightByReadsNoLowSubs	dnMinusDsWeightBySubsNoLowSubs
scatterplot_actual_v_expected("Expected", "dN/dS", "dNdSWeightByReads", "Site Ave dN/dS Weighted by Reads")
scatterplot_actual_v_expected("Expected", "dN/dS", "dNdSWeightByReadsNoLowSub", "Site Ave dN/dS Weighted by Reads \n(Exclude Windows with Low Site Substitutions)")
scatterplot_actual_v_expected("Expected", "dN/dS", "dNdSWeightBySubs", "Site Ave dN/dS Weighted by Substitutions")
scatterplot_actual_v_expected("Expected", "dN/dS", "dNdSWeightBySubsNoLowSub", "Site Ave dN/dS Weighted by Substitutions \n(Exclude Windows with Low Site Substitutions)")
scatterplot_actual_v_expected("ExpectedMinus", "dN-dS/TreeLen", "dnMinusDsWeightByReads", "Site Ave dN-dS/TreeLen Weighted By Reads")
scatterplot_actual_v_expected("ExpectedMinus", "dN-dS/TreeLen", "dnMinusDsWeightByReadsNoLowSubs", "Site Ave dN-dS/TreeLen Weighted By Reads\n(Exclude Windows with Low Site Substitutions)")
scatterplot_actual_v_expected("ExpectedMinus", "dN-dS/TreeLen", "dnMinusDsWeightBySubs", "Site Ave dN-dS/TreeLen Weighted By Substitutions")
scatterplot_actual_v_expected("ExpectedMinus", "dN-dS/TreeLen", "dnMinusDsWeightBySubsNoLowSubs", "Site Ave dN-dS/TreeLen Weighted By Substitutions\n(Exclude Windows with Low Site Substitutions)")

#' **Smoothed Line Plot Weighting Mechanisms By Site**
#' 
#+ fig.width=28
fullDatBydatsource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars=c("Ref", "Site"),
                                              measure.vars=c("dNdSWeightByReads", "dNdSWeightByReadsNoLowSub", 
                                                             "dNdSWeightBySubs", "dNdSWeightBySubsNoLowSub",
                                                             "dnMinusDsWeightByReads", "dnMinusDsWeightByReadsNoLowSubs",
                                                             "dnMinusDsWeightBySubs", "dnMinusDsWeightBySubsNoLowSubs",
                                                             "Expected", "ExpectedMinus"),
                                              variable.name="datsource", value.name="dnds")
head(fullDatBydatsource)
tail(fullDatBydatsource)
str(fullDatBydatsource)
summary(fullDatBydatsource)
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReads", "dNdSWeightBySubs", 
                                                              "dNdSWeightByReadsNoLowSub", "dNdSWeightBySubsNoLowSub", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReadsNoLowSub", "dNdSWeightBySubsNoLowSub", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Average Site dn/ds Weighted by Reads, Substitutions across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReads", "dNdSWeightBySubs", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Line Plot of Average Site dn/ds Weighted by Reads, Substitutions (Exclude Window Sites with Low Substitutions)**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReadsNoLowSub", "dNdSWeightBySubsNoLowSub", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Smoothed Line Plot of Ave Site dn-ds/Treelen across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dnMinusDsWeightByReads", "dnMinusDsWeightByReadsNoLowSubs", 
                                                              "dnMinusDsWeightBySubs", "dnMinusDsWeightBySubsNoLowSubs", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN-dS") + 
  ggtitle("Smoothed dn-ds") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Smoothed Line Plot of Ave Site dn-ds/Treelen across the genome, (Exclude Window Sites with Low Substitutions)**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dnMinusDsWeightByReadsNoLowSubs", 
                                                              "dnMinusDsWeightBySubsNoLowSubs", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN-dS/TreeLen") + 
  ggtitle("Smoothed dN-dS/Treelen") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Ave Site dn-ds/Treelen Weighted by Reads, Substitutions**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dnMinusDsWeightByReads", "dnMinusDsWeightBySubs", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN-dS/TreeLen") + 
  ggtitle("dN-dS/TreeLen by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Ave Site dn-ds/Treelen (Exclude Window Sites with Low Substitutions)**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dnMinusDsWeightByReadsNoLowSubs", "dnMinusDsWeightBySubsNoLowSubs", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN-dS") + 
  ggtitle("dN-dS/Treelen by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Plot the Ave Depth of Unambiguous Codons By Site**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=Codons) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Sequences with Unambiguous Codons") + 
  ggtitle("Codon Coverage Across Genome")

#' **Plot the Ave nonsynonymous substitutions across phylogeny**
#' 
fullDat$ExpectedSubs <- fullDat$ExpectedNonSyn + fullDat$ExpectedSyn
fullDatSubs <- reshape2:::melt(data=fullDat[, c("Site", "ExpectedNonSyn", "ExpectedSyn", "ExpectedSubs", "NonSyn", "Syn", "Subs")], 
                                id.vars="Site",
                                measure.vars=c("ExpectedNonSyn", "ExpectedSyn", "ExpectedSubs", "NonSyn", "Syn", "Subs"),
                                variable.name="subsSource", value.name="total")
summary(fullDatSubs)
head(fullDatSubs)

#+ fig.width=20
ggplot(fullDatSubs[fullDatSubs$subsSource %in% c("NonSyn", "ExpectedNonSyn"),], aes(x=Site, y=total, color=subsSource) ) + 
  geom_line() + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("Nonsynonymous Substitutions") + 
  ggtitle("Nonsynonymous Substitutions Across Phylogeny")

#' **Plot the synonymous substitutions across phylogeny**
#' 
#+ fig.width=20
ggplot(fullDatSubs[fullDatSubs$subsSource %in% c("Syn", "ExpectedSyn"),], aes(x=Site, y=total, color=subsSource) ) + 
  geom_line() + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("Synonymous Substitutions") + 
  ggtitle("Synonymous Substitutions Across Phylogeny")

#' **Plot the substitutions across phylogeny**
#' 
#+ fig.width=20
ggplot(fullDatSubs[fullDatSubs$subsSource %in% c("Subs", "ExpectedSubs"),], aes(x=Site, y=total, color=subsSource) ) + 
  geom_line() + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("Substitutions") + 
  ggtitle("Substitutions Across Phylogeny")

#' **Plot the Windows across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=Windows) ) + 
  geom_line() + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site Along Genome") + 
  ylab("Windows") + 
  ggtitle("Windows Across Genome")

#' TODO:  Strip indelible comparisons to another R file
#' 
#' **Plot the expected mutation scaling factor across the genome**
#' 
#+ fig.width=20
ggplot(indelible_dnds, aes(x=Site, y=as.numeric(as.character(Scaling_factor))) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Mutation Rate Scaling Factor") + 
  ggtitle("Mutation Scaling Factor Along Genome")

#' **Plot the Expected Omega rate across the genome**
#' 
#+ fig.width=20
ggplot(expected_dnds, aes(x=Site, y=Omega) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dn/dS Expected") + 
  ggtitle("Expected Selection Along Genome")


#' Concordance
#' ======================================
#' 

# Returns a table of Lin's concordance correlation values
print_table_corr <- function() {
  site_dnds_corr <- sapply(c("dNdSWeightByReads", "dNdSWeightBySubs", "dNdSWeightByReadsNoLowSub", "dNdSWeightBySubsNoLowSub"),
                           function(col) {
                             dnds_ccc <- epi.ccc(fullDat[, col], fullDat$Expected)
                             return (dnds_ccc$rho.c$est)
                           })

  
  site_dn_minus_ds_corr <- sapply(c("dnMinusDsWeightByReads", "dnMinusDsWeightBySubs", "dnMinusDsWeightByReadsNoLowSubs", "dnMinusDsWeightBySubsNoLowSubs"),
                                  function(col) {
                                    dnds_ccc <- epi.ccc(fullDat[, col], fullDat$ExpectedMinus)
                                    return (dnds_ccc$rho.c$est)
                                  })
  
  corr_vals <- data.frame(Concordance=c(site_dnds_corr,  site_dn_minus_ds_corr))
  corr_vals$Metric <- rownames(corr_vals)
  return (corr_vals)
}

table_corr <- print_table_corr()
kable(table_corr, format="html", caption="Concordance Correlation")

# Write out the concordance for validating in python unit tests
write.table(table_corr, file="./umberjack_unit_test.concordance.csv", row.names=FALSE, sep=",")

kable(table_corr, format="html", caption="Concordance Correlation")

#' **Find concordance at each mutation rate**
#' 
get_mut_concord <- function(scaling_factor) {
  ccc_dNdSWeightByReadsNoLowSub <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dNdSWeightByReadsNoLowSub, 
                       fullDat[fullDat$Scaling_factor==scaling_factor, ]$Expected)
  ccc_dNdSWeightByReads <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dNdSWeightByReads, 
                                            fullDat[fullDat$Scaling_factor==scaling_factor, ]$Expected)
  ccc_dNdSWeightBySubs <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dNdSWeightBySubs, 
                                    fullDat[fullDat$Scaling_factor==scaling_factor, ]$Expected)
  ccc_dNdSWeightBySubsNoLowSub <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dNdSWeightBySubsNoLowSub, 
                                   fullDat[fullDat$Scaling_factor==scaling_factor, ]$Expected)
  
  ccc_dnMinusDsWeightByReads <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dnMinusDsWeightByReads, 
                                                  fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedMinus)
  ccc_dnMinusDsWeightByReadsNoLowSubs <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dnMinusDsWeightByReadsNoLowSubs, 
                              fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedMinus)
  ccc_dnMinusDsWeightBySubs <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dnMinusDsWeightBySubs, 
                                         fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedMinus)
  ccc_dnMinusDsWeightBySubsNoLowSubs <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dnMinusDsWeightBySubsNoLowSubs, 
                                        fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedMinus)
  return (data.frame(Scaling_factor=scaling_factor, 
                     ccc_dNdSWeightByReadsNoLowSub$rho.c$est, 
                     ccc_dNdSWeightByReads$rho.c$est,
                     ccc_dNdSWeightBySubs$rho.c$est,
                     ccc_dNdSWeightBySubsNoLowSub$rho.c$est,
                     ccc_dnMinusDsWeightByReads$rho.c$est,
                     ccc_dnMinusDsWeightByReadsNoLowSubs$rho.c$est,
                     ccc_dnMinusDsWeightBySubs$rho.c$est,
                     ccc_dnMinusDsWeightBySubsNoLowSubs$rho.c$est))
}
dnds_ccc <- adply(.data=levels(fullDat$Scaling_factor), .margins=1, .fun=get_mut_concord)
dnds_ccc <- dnds_ccc[, -1]  # remove extraneous column with rownames

#+ results='asis'
kable(dnds_ccc, format="html", row.names=FALSE, caption="Concordance Correlation At Each Mutation Rate")


#' **Plot Correlation at each mutation rate**
#' 
plot_mut_corr<- function(scaling_factor) {
  fig <- ggplot(fullDat,
         aes(x=Expected, y=dNdSWeightByReadsNoLowSub)) + 
    geom_point(alpha=0.5, na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Scale ", scaling_factor, ": Expected dN/dS vs Inferred dn/ds Weighted by Reads\n(Excl Window Sites Low Subs) "))
  print(fig)
  
  fig <- ggplot(fullDat, 
                aes(x=Expected, y=dnMinusDsWeightByReadsNoLowSubs)) + 
    geom_point(alpha=0.5, na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Scale ", scaling_factor, ": Expected dN-dS/Tree vs Inferred dn-ds/Tree Weighted by Reads\n(Excl Window Sites Low Subs"))
  print(fig)
}

figs <- sapply(unique(fullDat$Scaling_factor), plot_mut_corr)


#' **Find substitutions at each mutation rate**

get_mut_subs <- function(scaling_factor) {
  return (data.frame(Scaling_factor=scaling_factor, 
                     AveExpectedSubs=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedSubs, na.rm=TRUE),
                     AveExpectedSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedSyn, na.rm=TRUE),
                     AveExpectedNonSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedNonSyn, na.rm=TRUE),
                     AveSubs=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Subs, na.rm=TRUE),
                     AveNonSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$NonSyn, na.rm=TRUE),
                     AveSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Syn, na.rm=TRUE),
                     AveCodons=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Codons, na.rm=TRUE),
                     AveWindows=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Windows, na.rm=TRUE)
                     
          ))
}
mut_subs <- adply(.data=levels(fullDat$Scaling_factor), .margins=1, .fun=get_mut_subs)
mut_subs <- mut_subs[, -1]  # remove extraneous column with rownames

#+ results='asis'
kable(mut_subs, format="html", row.names=FALSE, caption="Substitutions At Each Mutation Rate")

#' **Boxplot Substitutions at each Indelible Mutation Scaling Rate**
#' 
fullDatIndelibleCmpBySource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars=c("Ref", "Site", "Scaling_factor"),
                                                 measure.vars=c("Codons", "Windows", "NonSyn", "Syn",
                                                                "Subs", "ExpectedSubs", 
                                                                "ExpectedSyn", "ExpectedNonSyn"),
                                                 variable.name="datsource", value.name="val")
head(fullDatIndelibleCmpBySource)
summary(fullDatIndelibleCmpBySource)
#+ fig.width=12
fig <- ggplot(fullDatIndelibleCmpBySource, aes(x=Scaling_factor, y=val, color=datsource)) +
  geom_boxplot() + 
  xlab("\n Indelible Mutation Scaling Factor") + 
  ylab("Substitutions and Coverage Related Counts\n") + 
  ggtitle("Boxplot Substitutions and Counts By Indelible Mutation Scaling Factor")
print(fig)