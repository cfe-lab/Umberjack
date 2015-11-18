library(knitr)
opts_chunk$set(progress=FALSE, verbose=FALSE, warning=FALSE, message=FALSE, width=1500)
options(width=100)

library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)
library(plyr)

source("plot_helper.R")

# Read locations of input files from local umberjack_unit_test.config file
CONFIG_FILENAME <- "./umberjack_unit_test.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

actual_dnds_filename <- config[config$key=="ACTUAL_DNDS_FILENAME",]$val
expected_dnds_filename <-  config[config$key=="EXPECTED_DNDS_FILENAME",]$val


#'  ACTUAL_DNDS_FILENAME=`r actual_dnds_filename`
#'  
#'  EXPECTED_DNDS_FILENAME=`r expected_dnds_filename`
#'  

actual_dnds <- read.table(actual_dnds_filename, header=TRUE, na.strings="None", comment.char = "#", sep=",")
dim(actual_dnds)
head(actual_dnds)
str(actual_dnds)
summary(actual_dnds)

# Cols: Observed S Changes  Observed NS Changes	E[S Sites]	E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS	dn/ds
# Parse the expected dnds filename for it start and end nucleotide positions (1-based)

expected_dnds <- read.table(expected_dnds_filename, header=TRUE, sep="\t")  # Site  Interval	ScalingRate	Rate_class	Omega
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


#' Compare Codon Site dN/dS Averaged Over Windows With Full Population HyPhy dN/dS
#' ============================================================================
#' 
#' **Scatterplot actual vs expected dn ds together**

fullDat <- merge(x=actual_dnds, 
                 y=expected_dnds[, c("Site", "Omega", "Scaled.dN.dS", "Observed.S.Changes", "Observed.NS.Changes")], 
                 by="Site", all=FALSE, sort=TRUE)
colnames(fullDat)[grep("Omega", colnames(fullDat))] <- "Expected"
colnames(fullDat)[grep("Scaled.dN.dS", colnames(fullDat))] <- "ExpectedMinus"
colnames(fullDat)[grep("Observed.S.Changes", colnames(fullDat))] <- "ExpectedSyn"
colnames(fullDat)[grep("Observed.NS.Changes", colnames(fullDat))] <- "ExpectedNonSyn"
summary(fullDat)

#' **Scatterplot of Site dn/ds across the genome**
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
#+ fig.width=20
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

fullDat$ExpectedSubs <- fullDat$ExpectedNonSyn + fullDat$ExpectedSyn
fullDatSubs <- reshape2:::melt(data=fullDat[, c("Site", "ExpectedNonSyn", "ExpectedSyn", "ExpectedSubs", "NonSyn", "Syn", "Subs")], 
                               id.vars="Site",
                               measure.vars=c("ExpectedNonSyn", "ExpectedSyn", "ExpectedSubs", "NonSyn", "Syn", "Subs"),
                               variable.name="subsSource", value.name="total")
summary(fullDatSubs)
head(fullDatSubs)

#' **What are the general patterns of dnds by site?**
#' 
#+ fig.width=17
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReads", "dNdSWeightBySubs", 
                                                              "dNdSWeightByReadsNoLowSub", "dNdSWeightBySubsNoLowSub", "Expected"),], 
  aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank(),
        legend.position="bottom")

#+ fig.width=17
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReadsNoLowSub", "dNdSWeightBySubsNoLowSub", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank(),
        legend.position="bottom")


#' **How does dN/dS by Site relate to Codon Coverage, Window Coverage, Substitutions by Site?**

#' Line Plot of Average Site dn/ds Weighted by Reads, Substitutions across the genome
#' 
#+ fig.width=20, fig.height=12
dndsfig <- ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReads", "dNdSWeightBySubs", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))

codonfig <- ggplot(actual_dnds, aes(x=Site, y=Codons) ) + geom_line() + 
  ylab("Ave Win Site\nUnambiguous Codon Depth") +   
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))

nsfig <- ggplot(fullDatSubs[fullDatSubs$subsSource %in% c("NonSyn", "ExpectedNonSyn"),], aes(x=Site, y=total, color=subsSource) ) + 
  geom_line() + 
  ylab("Ave Window Site\nNonsyn Subs") + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))

sfig <- ggplot(fullDatSubs[fullDatSubs$subsSource %in% c("Syn", "ExpectedSyn"),], aes(x=Site, y=total, color=subsSource) ) + 
  geom_line() + 
  ylab("Ave Window Site\nSyn Subs") + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))

winfig <- ggplot(actual_dnds, aes(x=Site, y=Windows) ) + 
  geom_line() + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  ylab("Windows") + 
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))
  
list_gps <- AlignPlots(dndsfig, codonfig, nsfig, sfig, winfig)
do.call(grid.arrange, args=c(list_gps, ncol=1))

#' **Line Plot of Average Site dn/ds Weighted by Reads, Substitutions (Exclude Window Sites with Only Ambiguous Substitutions)**
#' 
#' Ambiguous substitutions refers to a site where there are only substitutions between ambiguous codons at the tip sequences.
#' 
dndsfig <- ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReadsNoLowSub", "dNdSWeightBySubsNoLowSub", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() +   
  ylab("dN/dS \nExclude Window Sites With Only Ambig Subs") + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))

#+ fig.width=20, fig.height=12
list_gps <- AlignPlots(dndsfig, codonfig, nsfig, sfig, winfig)
do.call(grid.arrange, args=c(list_gps, ncol=1))


#' **What are general patterns of dn-ds by site?**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dnMinusDsWeightByReads", "dnMinusDsWeightByReadsNoLowSubs", 
                                                              "dnMinusDsWeightBySubs", "dnMinusDsWeightBySubsNoLowSubs", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dn-ds/Treelen") + 
  ggtitle("Smoothed dn-ds/Treelen") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank(),
        legend.position="bottom")

#' **Smoothed Line Plot of Ave Site dn-ds/Treelen across the genome, (Exclude Window Sites with Only  Ambiguous Substitutions)**
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
        legend.text=element_text(size=24), legend.title=element_blank(),
        legend.position="bottom")


#' **How does dn-ds relate to Codons, Substitutions, Windows by site?**
#' 
dndsfig <- ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dnMinusDsWeightByReads", "dnMinusDsWeightBySubs", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  ylab("dN-dS/TreeLen") + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))

#+ fig.width=20, fig.height=12
list_gps <- AlignPlots(dndsfig, codonfig, nsfig, sfig, winfig)
do.call(grid.arrange, args=c(list_gps, ncol=1))

#' **Line Plot of Ave Site dn-ds/Treelen (Exclude Window Sites with Low Substitutions)**
#' 
dndsfig <- ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dnMinusDsWeightByReadsNoLowSubs", "dnMinusDsWeightBySubsNoLowSubs", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  ylab("dN-dS") + 
  scale_x_continuous(breaks = seq(0,  max(fullDatBydatsource$Site), 50)) + 
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"))

#+ fig.width=20, fig.height=12
list_gps <- AlignPlots(dndsfig, codonfig, nsfig, sfig, winfig)
do.call(grid.arrange, args=c(list_gps, ncol=1))

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

#' **Find concordance at varying substitution rates**
#' 
get_mut_concord <- function(x) {

  return (data.frame(dNdSWeightByReadsNoLowSub = epi.ccc(x$dNdSWeightByReadsNoLowSub, x$Expected)$rho.c$est, 
                     dNdSWeightByReads = epi.ccc(x$dNdSWeightByReads, x$Expected)$rho.c$est,
                     dNdSWeightBySubs = epi.ccc(x$dNdSWeightBySubs, x$Expected)$rho.c$est,
                     dNdSWeightBySubsNoLowSub =  epi.ccc(x$dNdSWeightBySubsNoLowSub, x$Expected)$rho.c$est,
                     dnMinusDsWeightByReads = epi.ccc(x$dnMinusDsWeightByReads, x$ExpectedMinus)$rho.c$est,
                     dnMinusDsWeightByReadsNoLowSubs = epi.ccc(x$dnMinusDsWeightByReadsNoLowSubs, x$ExpectedMinus)$rho.c$est,
                     dnMinusDsWeightBySubs = epi.ccc(x$dnMinusDsWeightBySubs,  x$ExpectedMinus)$rho.c$est,
                     dnMinusDsWeightBySubsNoLowSubs = epi.ccc(x$dnMinusDsWeightBySubsNoLowSubs, x$ExpectedMinus)$rho.c$est))
}

# Bin the sites by the number of site substitutions (counted across phylogeny)
fullDat$SubBin <- cut(fullDat$Subs, breaks=c(0, 1, 5, 10, 20, 50, 100))
dnds_ccc <- ddply(.data=fullDat, 
                  .variables="SubBin",
                  .fun=get_mut_concord)


#+ results='asis'
kable(dnds_ccc, format="html", row.names=FALSE, caption="Concordance Correlation At Each Mutation Rate")


#' **Plot Correlation at each mutation rate**
#' 
plot_mut_corr<- function(ScalingRate) {
  
  if (nrow(fullDat[!is.na(fullDat$SubBin) & fullDat$SubBin == ScalingRate,]) == 0) {
    return ()
  }
  
  fig <- ggplot(fullDat[fullDat$SubBin == ScalingRate,],
         aes(x=Expected, y=dNdSWeightByReadsNoLowSub)) + 
    geom_point(alpha=0.5, na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Scale ", ScalingRate, ": Expected dN/dS vs Inferred dn/ds Weighted by Reads\n(Excl Window Sites Low Subs) "))
  print(fig)
  
  fig <- ggplot(fullDat[fullDat$SubBin == ScalingRate,],
                aes(x=Expected, y=dnMinusDsWeightByReadsNoLowSubs)) + 
    geom_point(alpha=0.5, na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Scale ", ScalingRate, ": Expected dN-dS vs Inferred dn-ds Weighted by Reads\n(Excl Window Sites Low Subs"))
  print(fig)
}

figs <- sapply(levels(fullDat$SubBin), plot_mut_corr)

