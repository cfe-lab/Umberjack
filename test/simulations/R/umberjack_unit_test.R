library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)
library(plyr)

NUC_PER_CODON <- 3


# Read locations of input files from local sliding_window_tree_unit_test.config file
CONFIG_FILENAME <- "./sliding_window_tree_unit_test.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

actual_dnds_filename <- config[config$key=="ACTUAL_DNDS_FILENAME",]$val
expected_dnds_filename <-  config[config$key=="EXPECTED_DNDS_FILENAME",]$val



indelible_dnds_filename <- config[config$key=="INDELIBLE_DNDS_FILENAME", ]$val

#'  ACTUAL_DNDS_FILENAME=`r actual_dnds_filename`
#'  
#'  EXPECTED_DNDS_FILENAME=`r expected_dnds_filename`
#'  
#'  INDELIBLE_DNDS_FILENAME=`r indelible_dnds_filename`

actual_dnds_file <- file(actual_dnds_filename, open="rt")
comments <- readLines(actual_dnds_file, 1) # Read one line 
close(actual_dnds_file)

args <- unlist(strsplit(comments, ','))
smooth_dist <- as.numeric(unlist(strsplit(args[grep("smooth_dist", args)], "="))[2])

#'  smooth_dist=`r smooth_dist`

#' **`r comments`**
#' -----------------------------
#' 
#' 
# Cols:  Ref  Site  aveDnDs  dNdSWeightBySubst	dN_minus_dS	Windows	Codons	NonSyn	Syn	Subst	dNdSWeightByReads	multisiteAvedNdS	multisitedNdSWeightBySubst  dNdSWeightByReadsNoLowSyn
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
dim(expected_dnds)
head(expected_dnds)
str(expected_dnds)
summary(expected_dnds)


# check consistency
all(expected_dnds$Site == actual_dnds$Site)


MIN_SITE <- min(expected_dnds$Site)
MAX_SITE <- max(expected_dnds$Site)
expected_dnds$MultisiteAveDnDs <- apply(expected_dnds, 1, function(row) {
  site <- as.numeric(row["Site"])
  total_dnds <- 0
  total_sites <- 0
  for (adjsite in max(MIN_SITE, site-smooth_dist):min(MAX_SITE, site+smooth_dist)) {
    if (!is.na(expected_dnds[expected_dnds$Site==adjsite, "Omega"])) {
      total_dnds <- total_dnds + expected_dnds[expected_dnds$Site==adjsite, "Omega"]
      total_sites <- total_sites +1
    }
  }
  multisite_ave <- total_dnds / total_sites
  return (multisite_ave)
})


#' Compare Full Population INDELible dN/dS with Full Population HyPhy dN/dS
#' ============================================================================

# Read in Indelible's intended dN/dS
# Cols: Site,Interval,Scaling_factor,Rate_class,Omega
indelible_dnds <- read.table(indelible_dnds_filename, header=TRUE, sep=",")
summary(indelible_dnds)
head(indelible_dnds)

#' **Plot indelible's intended dN/dS vs the actual dN/dS as determined by HyPhy**
#' 
indelhyphy <- data.frame(Site=indelible_dnds$Site, 
                    Scaling_factor=indelible_dnds$Scaling_factor, 
                    Rate_class=indelible_dnds$Rate_class,
                    Indelible_Omega=indelible_dnds$Omega,
                    Hyphy_Omega=expected_dnds$Omega)
summary(indelhyphy)
indelhyphy <- na.omit(indelhyphy)

#' **Scatterplot Full Population HyPhy Actual dN/dS vs Indelible Expected dN/dS**
ggplot(indelhyphy, aes(x=Indelible_Omega, y=Hyphy_Omega)) + 
  geom_point(shape=1, alpha=0.5, na.rm=TRUE) +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
  geom_abline(slope=1) + 
  ylab("HyPhy dN/dS\n") + 
  xlab("\nIndelible dN/dS") + 
  #coord_fixed(ratio=1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  ggtitle("HyPhy dN/dS Vs Indelible dN/dS on Full Population")

#' **Boxplot Full Population HyPhy Actual dN/dS vs Indelible Expected dN/dS**
ggplot(indelhyphy, aes(x=as.factor(Indelible_Omega), y=Hyphy_Omega)) + 
  geom_boxplot() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
  ylab("HyPhy dN/dS\n") + 
  xlab("\nIndelible dN/dS") + 
  scale_y_continuous(breaks=seq(0, max(indelhyphy$Hyphy_Omega), 0.5)) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=14) ,
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("HyPhy dN/dS Vs Indelible dN/dS on Full Population")



#' Compare Codon Site dN/dS Averaged Over Windows With Full Population HyPhy dN/dS
#' ============================================================================
#' 
#' **Scatterplot actual vs expected dn ds together**
fullDat <- cbind(actual_dnds, 
                 Expected=expected_dnds$Omega, 
                 ExpectedMultisite=expected_dnds$MultisiteAveDnDs, 
                 ExpectedMinus=expected_dnds$Scaled.dN.dS,
                 ExpectedSyn=expected_dnds$Observed.S.Changes,
                 ExpectedNonSyn=expected_dnds$Observed.NS.Changes,
                 Scaling_factor=indelible_dnds$Scaling_factor)
summary(fullDat)


scatterplot_actual_v_expected <- function(expected_colname, expected_title, actual_colname, actual_title) {
  fig <- ggplot(fullDat  , aes_string(x=expected_colname, y=actual_colname)) + 
    geom_point(na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Inferred ", actual_title, " Vs Expected ", expected_colname) )
  print(fig)
  
  actual_col <- fullDat[, actual_colname]
  ave <- mean(actual_col, na.rm=TRUE) 
  most_outlier <- max(abs( c(max(actual_col - ave, na.rm=TRUE), 
                             min(actual_col - ave, na.rm=TRUE))))
  stdev <- sd(actual_col, na.rm=TRUE)
  if (most_outlier > 4 * stdev) {
    fig <- ggplot(fullDat[abs(actual_col - ave) < 3 * stdev, ]  , aes_string(x=expected_colname, y=actual_colname)) + 
      geom_point(na.rm=TRUE) +
      geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
      geom_abline(slope=1) + 
      ylab("Inferred\n") + 
      xlab("\nExpected") + 
      #coord_fixed(ratio=1) + 
      theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
      ggtitle(paste0("No Outliers: Inferred ", actual_title, " Vs Expected ", expected_title) )
    print(fig)
  }
}

scatterplot_actual_v_expected("Expected", "dN/dS", "aveDnDs", "Site Ave dN/dS")
scatterplot_actual_v_expected("Expected", "dN/dS", "dNdSWeightByReads", "Site Ave dN/dS Weighted by Reads")
scatterplot_actual_v_expected("Expected", "dN/dS", "dNdSWeightByReadsNoLowSyn", "Site Ave dN/dS Weighted by Reads (Exclude Low Syn Subst)")
scatterplot_actual_v_expected("Expected", "dN/dS", "dNdSWeightBySubst", "Site Ave dN/dS Weighted by Substitutions")
scatterplot_actual_v_expected("ExpectedMinus", "dN-dS", "dN_minus_dS", "Site Ave dN-dS")
scatterplot_actual_v_expected("ExpectedMultisite", "Multisite Ave dN/dS", "multisiteAvedNdS", "Multisite Ave dN/dS")
scatterplot_actual_v_expected("ExpectedMultisite", "Multisite Ave dN/dS", "multisitedNdSWeightBySubst", "Multisite Ave dN/dS Weighted by Substitutions")


#' **Smoothed Scatterplot of Site dn/ds across the genome**
#+ fig.width=28
fullDatBydatsource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars=c("Ref", "Site"),
                                              measure.vars=c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "dN_minus_dS",
                                                             "multisiteAvedNdS", "multisitedNdSWeightBySubst", 
                                                             "Expected", "ExpectedMultisite", "ExpectedMinus", 
                                                             "dNdSWeightByReadsNoLowSyn"),
                                              variable.name="datsource", value.name="dnds")
head(fullDatBydatsource)
tail(fullDatBydatsource)
str(fullDatBydatsource)
summary(fullDatBydatsource)
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReadsNoLowSyn", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Average Site dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("aveDnDs", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Line Plot of Average Site dn/ds Weighted by Substitutions and Reads across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightBySubst", "dNdSWeightByReads", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Line Plot of Average Site dn/ds Weighted by Reads With Windows with < 1 Syn Subst removed across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dNdSWeightByReadsNoLowSyn", "Expected"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Site dn-ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dN_minus_dS", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN-dS") + 
  ggtitle("dn-ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Smoothed Scatterplot of Site dn-ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("dN_minus_dS", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN-dS") + 
  ggtitle("Smoothed dn-ds") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Smoothed Scatterplot of MultiSite dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("multisiteAvedNdS", "multisitedNdSWeightBySubst", "ExpectedMultisite"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("Multisite dn/ds") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of MultiSite dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBydatsource[fullDatBydatsource$datsource %in% c("multisiteAvedNdS", "multisitedNdSWeightBySubst", "ExpectedMultisite"),], 
       aes(x=Site, y=dnds, color=datsource) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("Multisite dn/ds") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())



#' **Plot the Actual Codon Coverage across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=Codons) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Sequences with Unambiguous Codons") + 
  ggtitle("Codon Coverage Across Genome")

#' **Plot the nonsynonymous substitutions across phylogeny**
#' 
fullDat$ExpectedSubst <- fullDat$ExpectedNonSyn + fullDat$ExpectedSyn
fullDatSubst <- reshape2:::melt(data=fullDat[, c("Site", "ExpectedNonSyn", "ExpectedSyn", "ExpectedSubst", "NonSyn", "Syn", "Subst")], 
                                id.vars="Site",
                                measure.vars=c("ExpectedNonSyn", "ExpectedSyn", "ExpectedSubst", "NonSyn", "Syn", "Subst"),
                                variable.name="substSource", value.name="total")
summary(fullDatSubst)
head(fullDatSubst)

#+ fig.width=20
ggplot(fullDatSubst[fullDatSubst$substSource %in% c("NonSyn", "ExpectedNonSyn"),], aes(x=Site, y=total, color=substSource) ) + 
  geom_line() + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("Nonsynonymous Substitutions") + 
  ggtitle("Nonsynonymous Substitutions Across Phylogeny")

#' **Plot the synonymous substitutions across phylogeny**
#' 
#+ fig.width=20
ggplot(fullDatSubst[fullDatSubst$substSource %in% c("Syn", "ExpectedSyn"),], aes(x=Site, y=total, color=substSource) ) + 
  geom_line() + 
  geom_smooth(na.rm=TRUE) + 
  xlab("Codon Site") + 
  ylab("Synonymous Substitutions") + 
  ggtitle("Synonymous Substitutions Across Phylogeny")

#' **Plot the substitutions across phylogeny**
#' 
#+ fig.width=20
ggplot(fullDatSubst[fullDatSubst$substSource %in% c("Subst", "ExpectedSubst"),], aes(x=Site, y=total, color=substSource) ) + 
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

#' TODO:  Strip indelible compjarisons to another R file
#' 
#' **Plot the expected mutation scaling factor across the genome**
#' 
#+ fig.width=20
ggplot(indelible_dnds, aes(x=Site, y=Scaling_factor) ) + geom_line() + 
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


# Returns a table of Lin's concordance correlation values
print_table_corr <- function() {
  site_dnds_corr <- sapply(c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "dNdSWeightByReadsNoLowSyn"),
                           function(col) {
                             dnds_ccc <- epi.ccc(fullDat[, col], fullDat$Expected)
                             return (dnds_ccc$rho.c$est)
                           })
  
  multisite_dnds_corr <- sapply(c("multisiteAvedNdS", "multisitedNdSWeightBySubst"),
                                function(col) {
                                  dnds_ccc <- epi.ccc(fullDat[, col], fullDat$ExpectedMultisite)
                                  return (dnds_ccc$rho.c$est)
                                })
  
  site_dn_minus_ds_corr <- sapply(c("dN_minus_dS"),
                                  function(col) {
                                    dnds_ccc <- epi.ccc(fullDat[, col], fullDat$ExpectedMinus)
                                    return (dnds_ccc$rho.c$est)
                                  })
  corr_vals <- data.frame(Concordance=c(site_dnds_corr, multisite_dnds_corr, site_dn_minus_ds_corr))
}

table_corr <- print_table_corr()
kable(table_corr, format="html", caption="Concordance Correlation")

# **Find concordance at each mutation rate**

get_mut_concord <- function(scaling_factor) {
  ccc_dNdSWeightByReadsNoLowSyn <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dNdSWeightByReadsNoLowSyn, 
                       fullDat[fullDat$Scaling_factor==scaling_factor, ]$Expected)
  ccc_dNdSWeightByReads <-  epi.ccc(fullDat[fullDat$Scaling_factor==scaling_factor, ]$dNdSWeightByReads, 
                                            fullDat[fullDat$Scaling_factor==scaling_factor, ]$Expected)
  return (data.frame(Scaling_factor=scaling_factor, ccc_dNdSWeightByReadsNoLowSyn$rho.c$est, ccc_dNdSWeightByReads$rho.c$est))
}
dnds_ccc <- adply(.data=unique(fullDat$Scaling_factor), .margins=1, .fun=get_mut_concord)
dnds_ccc <- dnds_ccc[, -1]  # remove extraneous column with rownames

#+ results='asis'
kable(dnds_ccc, format="html", row.names=FALSE, caption="Concordance Correlation At Each Mutation Rate")


# **Plot Correlation at each mutation rate**

plot_mut_corr<- function(scaling_factor) {
  fig <- ggplot(fullDat[fullDat$Scaling_factor==scaling_factor & fullDat$ExpectedSyn >=1, ],
         aes(x=Expected, y=dNdSWeightByReadsNoLowSyn)) + 
    geom_point(na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Scale ", scaling_factor, ": Expected dN/dS vs Inferred dn/ds Weighted by Subst (Exclude Low Subst) "))
  print(fig)
  
  fig <- ggplot(fullDat[fullDat$Scaling_factor==scaling_factor & fullDat$ExpectedSyn >=1, ], 
                aes(x=Expected, y=dNdSWeightByReads)) + 
    geom_point(na.rm=TRUE) +
    geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC", na.rm=TRUE) +
    geom_abline(slope=1) + 
    ylab("Inferred\n") + 
    xlab("\nExpected") + 
    #coord_fixed(ratio=1) + 
    theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Scale ", scaling_factor, ": Expected dN/dS vs Inferred dn/ds Weighted by Subst"))
  print(fig)
}

sapply(unique(fullDat$Scaling_factor), plot_mut_corr)


# **Find substitutions at each mutation rate**

get_mut_subst <- function(scaling_factor) {
  return (data.frame(Scaling_factor=scaling_factor, 
                     AveExpectedSubst=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedSubst, na.rm=TRUE),
                     AveExpectedSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedSyn, na.rm=TRUE),
                     AveExpectedNonSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$ExpectedNonSyn, na.rm=TRUE),
                     AveSubst=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Subst, na.rm=TRUE),
                     AveNonSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$NonSyn, na.rm=TRUE),
                     AveSyn=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Syn, na.rm=TRUE),
                     AveCodons=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Codons, na.rm=TRUE),
                     AveWindows=mean(fullDat[fullDat$Scaling_factor==scaling_factor, ]$Windows, na.rm=TRUE)
                     
          ))
}
mut_subst <- adply(.data=unique(fullDat$Scaling_factor), .margins=1, .fun=get_mut_subst)
mut_subst <- mut_subst[, -1]  # remove extraneous column with rownames

#+ results='asis'
kable(mut_subst, format="html", row.names=FALSE, caption="Substitutions At Each Mutation Rate")