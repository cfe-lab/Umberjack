library(limma)
library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)

BASE_MUTATION_RATE <- 0.0001  # TODO:  output mutation rate into comments so that we can parse it in R

NUC_PER_CODON <- 3


# Read locations of input files from local sliding_window_tree_unit_test.config file
CONFIG_FILENAME <- "./sliding_window_tree_unit_test.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

actual_dnds_filename <- config[config$key=="ACTUAL_DNDS_FILENAME",]$val
expected_dnds_filename <-  config[config$key=="EXPECTED_DNDS_FILENAME",]$val
expected_dnds_start_nuc_pos <-  as.numeric(config[config$key=="EXPECTED_DNDS_START_NUC_POS",]$val)
expected_dnds_end_nuc_pos <-  as.numeric(config[config$key=="EXPECTED_DNDS_END_NUC_POS",]$val)

#'  ACTUAL_DNDS_FILENAME=`r actual_dnds_filename`
#'  
#'  EXPECTED_DNDS_FILENAME=`r expected_dnds_filename`
#'  
#'  EXPECTED_DNDS_START_NUC_POS=`r expected_dnds_start_nuc_pos`
#'  
#'  EXPECTED_DNDS_END_NUC_POS=`r expected_dnds_end_nuc_pos`

dnds_file <- file(actual_dnds_filename, open="rt")
comments <- readLines(dnds_file, 1) # Read one line 
close(dnds_file)

args <- unlist(strsplit(comments, ','))
start_nuc_pos <- as.numeric(unlist(strsplit(args[grep("start_nuc_pos", args)], "="))[2])
end_nuc_pos <- as.numeric(unlist(strsplit(args[grep("end_nuc_pos", args)], "="))[2])
start_codon <- (start_nuc_pos %/% 3) + 1
end_codon <- end_nuc_pos %/% 3
smooth_dist <- as.numeric(unlist(strsplit(args[grep("smooth_dist", args)], "="))[2])

#'  start_codon=`r start_codon`
#'  
#'  end_codon=`r end_codon`
#'  
#'  smooth_dist=`r smooth_dist`

#' **`r comments`**
#' -----------------------------
#' 
#' 
# Cols:  Ref  Site  aveDnDs  dNdSWeightBySubst	dN_minus_dS	Windows	Codons	NonSyn	Syn	Subst	dNdSWeightByReads	multisiteAvedNdS	multisitedNdSWeightBySubst	simpleDnDs  dNdSWeightByReadsNoLowSyn
actual_dnds <- read.table(actual_dnds_filename, header=TRUE, na.strings="None", comment.char = "#", sep=",")
dim(actual_dnds)
head(actual_dnds)
str(actual_dnds)
summary(actual_dnds)

# Cols: Observed S Changes  Observed NS Changes	E[S Sites]	E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS	dn/ds
# Parse the expected dnds filename for it start and end nucleotide positions (1-based)

expected_dnds <- read.table(expected_dnds_filename, header=TRUE, sep="\t")  # Site  Interval	Scaling_factor	Rate_class	Omega
expected_dnds_start_codon <- (expected_dnds_start_nuc_pos %/% 3) + 1
expected_dnds$Site <- as.numeric(rownames(expected_dnds)) + expected_dnds_start_codon -1
expected_dnds_end_codon <- max(expected_dnds$Site)
expected_dnds$Omega <- expected_dnds$dN/expected_dnds$dS
expected_dnds$Omega[expected_dnds$dS == 0] <- NA
dim(expected_dnds)
head(expected_dnds)
str(expected_dnds)
summary(expected_dnds)


intersect_start_codon <- max(expected_dnds_start_codon, start_codon)
intersect_end_codon <- min(expected_dnds_end_codon, end_codon)
actual_dnds <- actual_dnds[actual_dnds$Site >= intersect_start_codon & actual_dnds$Site <=intersect_end_codon,]
summary(actual_dnds)
expected_dnds <- expected_dnds[expected_dnds$Site >= intersect_start_codon & expected_dnds$Site <=intersect_end_codon,]
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
    total_dnds <- total_dnds + expected_dnds[expected_dnds$Site==adjsite, "Omega"]
    total_sites <- total_sites +1
  }
  multisite_ave <- total_dnds / total_sites
  return (multisite_ave)
})





#' **Scatterplot actual vs expected dn ds together**
fullDat <- cbind(actual_dnds, Expected=expected_dnds$Omega, ExpectedMultisite=expected_dnds$MultisiteAveDnDs, ExpectedMinus=expected_dnds$Scaled.dN.dS)
fullDat <- na.omit(fullDat)
summary(fullDat)

ggplot(fullDat[!is.na(fullDat$aveDnDs),], aes(x=Expected, y=aveDnDs)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  geom_abline(slope=1) + 
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  coord_fixed(ratio=1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  ggtitle("Inferred Ave dN/dS vs Expected dn/ds")


ggplot(fullDat[!is.na(fullDat$simpleDnDs),], aes(x=Expected, y=simpleDnDs)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  coord_fixed(ratio=1) + 
  geom_abline(slope=1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  ggtitle("Inferred Simple dN/dS (do not normalize by expected subst) vs Expected dn/ds")


ggplot(fullDat[!is.na(fullDat$dNdSWeightByReadsNoLowSyn),], aes(x=Expected, y=dNdSWeightByReadsNoLowSyn)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  coord_fixed(ratio=1) + 
  geom_abline(slope=1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  ggtitle("Inferred dN/dS weighted by reads (windows with low syn subst excl) vs Expected dn/ds")

ggplot(fullDat[!is.na(fullDat$dN_minus_dS),], aes(x=ExpectedMinus, y=dN_minus_dS)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  coord_fixed(ratio=1) + 
  geom_abline(slope=1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  ggtitle("Inferred dN-dS vs Expected dn-dS")

ggplot(fullDat[!is.na(fullDat$dNdSWeightBySubst),], aes(x=Expected, y=dNdSWeightBySubst)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  geom_abline(slope = 1) + 
  coord_fixed(ratio=1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  ggtitle("Inferred dn/ds Weighted by Substitutions vs Expected dn/ds")

ggplot(fullDat[!is.na(fullDat$dNdSWeightByReads),], aes(x=Expected, y=dNdSWeightByReads)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  geom_abline(slope = 1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  coord_fixed(ratio=1) + 
  ggtitle("Inferred dN/dS Weighted by Reads vs Expected dn/ds")

ggplot(fullDat, aes(x=ExpectedMultisite, y=multisiteAvedNdS)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  geom_abline(slope = 1) + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  coord_fixed(ratio=1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Averaged Across Surrounding Sites")

ggplot(fullDat, aes(x=ExpectedMultisite, y=multisitedNdSWeightBySubst)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  coord_fixed(ratio=1) + 
  geom_abline(slope = 1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Averaged Across Surrounding Sites, Weighted by Substitutions")


#' **Smoothed Scatterplot of Site dn/ds across the genome**
#+ fig.width=28
fullDatBySource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars=c("Ref", "Site"),
                                              measure.vars=c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "dN_minus_dS",
                                                             "multisiteAvedNdS", "multisitedNdSWeightBySubst", 
                                                             "Expected", "ExpectedMultisite", "ExpectedMinus", 
                                                             "simpleDnDs", "dNdSWeightByReadsNoLowSyn"),
                                              variable.name="source", value.name="dnds")
head(fullDatBySource)
tail(fullDatBySource)
str(fullDatBySource)
summary(fullDatBySource)
ggplot(fullDatBySource[fullDatBySource$source %in% c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "Expected"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

ggplot(fullDatBySource[fullDatBySource$source %in% c("dNdSWeightByReadsNoLowSyn", "Expected"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Average Site dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("aveDnDs", "Expected"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Line Plot of Average Site dn/ds Weighted by Substitutions and Reads across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("dNdSWeightBySubst", "dNdSWeightByReads", "Expected"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Line Plot of Average Site dn/ds Weighted by Reads With Windows with < 1 Syn Subst removed across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("dNdSWeightByReadsNoLowSyn", "Expected"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Site dn-ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("dN_minus_dS", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN-dS") + 
  ggtitle("dn-ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Smoothed Scatterplot of Site dn-ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("dN_minus_dS", "ExpectedMinus"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("dN-dS") + 
  ggtitle("Smoothed dn-ds") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())


#' **Smoothed Scatterplot of MultiSite dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("multisiteAvedNdS", "multisitedNdSWeightBySubst", "ExpectedMultisite"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("Multisite dn/ds") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of MultiSite dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("multisiteAvedNdS", "multisitedNdSWeightBySubst", "ExpectedMultisite"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("Multisite dn/ds") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
        legend.text=element_text(size=24), legend.title=element_blank())

#' **Point Plot of dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource, aes(x=Site, y=dnds, color=source) ) + geom_point() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site")


#' **Plot the Codon Coverage across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=Codons) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Sequences with Unambiguous Codons") + 
  ggtitle("Codon Coverage Across Genome")

#' **Plot the nonsynonymous substitutions across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=NonSyn) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Nonsynonymous Substitutions") + 
  ggtitle("Total Genome Nonsynonymous Substitutions Across Population")

#' **Plot the synonymous substitutions across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=Syn) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Synonymous Substitutions") + 
  ggtitle("Total Genome Synonymous Substitutions Across Population")

#' **Plot the substitutions across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=Subst) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Substitutions") + 
  ggtitle("Total Genome Substitutions Across Population")

#' **Plot the Windows across genome**
#' 
#+ fig.width=20
ggplot(actual_dnds, aes(x=Site, y=Windows) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Windows") + 
  ggtitle("Windows Across Genome")

# #' **Plot the expected mutation rate across the genome**
# #' 
# #+ fig.width=20
# ggplot(expected_dnds, aes(x=Site, y=Scaling_factor*BASE_MUTATION_RATE) ) + geom_line() + 
#   xlab("Codon Site Along Genome") + 
#   ylab("Mutation Rate Scaling Factor") + 
#   scale_y_continuous(breaks=seq(0, 0.25, 0.01), limits=c(0, 0.25)) +
#   ggtitle("Mutation Rate Along Genome")

#' **Plot the Expected Omega rate across the genome**
#' 
#+ fig.width=20
ggplot(expected_dnds, aes(x=Site, y=Omega) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dn/dS Expected") + 
  ggtitle("Expected Selection Along Genome")

#' **Paired test without assumption of normalcy**
htest <-  wilcox.test(fullDat$dNdSWeightBySubst, fullDat$Expected, paired=TRUE, exact=TRUE, na.action="na.exclude")
print (htest)


#' **Find how correlated the actual ave dn/ds  vs expected dn/ds are**
dnds_cor <- cor(fullDat$aveDnDs, fullDat$Expected, method="spearman", use="complete.obs")
print(dnds_cor)


#' **Find how correlated the actual simple  dn/ds  vs expected dn/ds are**
dnds_cor <- cor(fullDat$simpleDnDs, fullDat$Expected, method="spearman", use="complete.obs")
print(dnds_cor)


#' **Find how correlated the actual dn/ds weighted by substitutions vs expected dn/ds are**
dnds_cor <- cor(fullDat$dNdSWeightBySubst, fullDat$Expected, method="spearman", use="complete.obs")
print(dnds_cor)

#' **Find how correlated the actual dn-ds vs expected dn-ds are**
dnds_cor <- cor(fullDat$dN_minus_dS, fullDat$ExpectedMinus, method="spearman", use="complete.obs")
print(dnds_cor)

#' **Find how correlated the actual dn/ds weighted by reads vs expected dn/ds are**
dnds_cor <- cor(fullDat$dNdSWeightByReads, fullDat$Expected, method="spearman", use="complete.obs")
print(dnds_cor)

#' **Find how correlated the actual dn/ds weighted by reads (windows low syn subs excl) vs expected dn/ds are**
dnds_cor <- cor(fullDat$dNdSWeightByReadsNoLowSyn, fullDat$Expected, method="spearman", use="complete.obs")
print(dnds_cor)


#' **Find how correlated the actual dn/ds smoothed over surrounding sites vs expected dn/ds are**
dnds_cor <- cor(fullDat$multisiteAvedNdS, fullDat$ExpectedMultisite, method="spearman", use="complete.obs")
print(dnds_cor)

# **Find how correlated the actual dn/ds smoothed over surrounding sites weighted by substitutions vs expected dn/ds are**
dnds_cor <- cor(fullDat$multisitedNdSWeightBySubst, fullDat$ExpectedMultisite, method="spearman", use="complete.obs")
print(dnds_cor)


#' **Find how concordance correlated the actual ave dn/ds  vs expected dn/ds are**
dnds_ccc <- epi.ccc(fullDat$aveDnDs, fullDat$Expected)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual simple  dn/ds  vs expected dn/ds are**
dnds_ccc <- epi.ccc(fullDat$simpleDnDs, fullDat$Expected)
print(dnds_ccc$rho.c)


#' **Find how concordance correlated the actual dn/ds weighted by substitutions vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$dNdSWeightBySubst, fullDat$Expected)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual dn-ds vs expected dn-ds are**
dnds_ccc <-  epi.ccc(fullDat$dN_minus_dS, fullDat$ExpectedMinus)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual dn/ds weighted by reads vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$dNdSWeightByReads, fullDat$Expected)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual dn/ds weighted by reads (windows low syn subst excl) vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$dNdSWeightByReadsNoLowSyn, fullDat$Expected)
print(dnds_ccc$rho.c)


#' **Find how concordance correlated the actual dn/ds smoothed over surrounding sites vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$multisiteAvedNdS, fullDat$ExpectedMultisite)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual dn/ds smoothed over surrounding sites weighted by substitutions vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$multisitedNdSWeightBySubst, fullDat$ExpectedMultisite)
print(dnds_ccc$rho.c)

# **Find concordance at each mtuation rate  5x
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site <= 400, ]$dNdSWeightByReadsNoLowSyn, fullDat[fullDat$Site <= 400, ]$Expected)
print(dnds_ccc$rho.c)

# 10x
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site >400 & fullDat$Site <=800 , ]$dNdSWeightByReadsNoLowSyn, fullDat[fullDat$Site > 400 & fullDat$Site <=800 , ]$Expected)
print(dnds_ccc$rho.c)

# 50x
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site >800 & fullDat$Site <=1200 , ]$dNdSWeightByReadsNoLowSyn, fullDat[fullDat$Site > 800 & fullDat$Site <=1200 , ]$Expected)
print(dnds_ccc$rho.c)

# 100x
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site >1200 & fullDat$Site <=1600 , ]$dNdSWeightByReadsNoLowSyn, fullDat[fullDat$Site > 1200 & fullDat$Site <=1600 , ]$Expected)
print(dnds_ccc$rho.c)

# **Find concordance at each mtuation rate  5x.  Don't ignore windows with low syn subs
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site <= 400, ]$dNdSWeightByReads, fullDat[fullDat$Site <= 400, ]$Expected)
print(dnds_ccc$rho.c)

# 10x
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site >400 & fullDat$Site <=800 , ]$dNdSWeightByReads, fullDat[fullDat$Site > 400 & fullDat$Site <=800 , ]$Expected)
print(dnds_ccc$rho.c)

# 50x
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site >800 & fullDat$Site <=1200 , ]$dNdSWeightByReads, fullDat[fullDat$Site > 800 & fullDat$Site <=1200 , ]$Expected)
print(dnds_ccc$rho.c)

# 100x
dnds_ccc <-  epi.ccc(fullDat[fullDat$Site >1200 & fullDat$Site <=1600 , ]$dNdSWeightByReads, fullDat[fullDat$Site > 1200 & fullDat$Site <=1600 , ]$Expected)
print(dnds_ccc$rho.c)
