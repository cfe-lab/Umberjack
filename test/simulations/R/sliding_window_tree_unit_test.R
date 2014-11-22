library(limma)
library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)

BASE_MUTATION_RATE <- 0.0001

NUC_PER_CODON <- 3
WINDOWSIZE <- 400
REF_LEN_NUC <- 9000

POPN_SIZE <- 10000
START_NUCPOS <- 3001
END_NUCPOS <- 3400

START_CODON <- 1001
END_CODON <- 1133

SMOOTH_DIST <- 10
# TODO:  find expected dn-ds


actual_dnds_filename <- "../out/consensus/mut100_1x_errfree_window400/actual_dnds_by_site.csv"
#actual_dnds_filename <- "../out/consensus/mut100_1x_window400/actual_dnds_by_site.csv"
expected_dnds_filename <- "../data/indelible/sample_genomes.100.rates.csv"

dnds_file <- file(actual_dnds_filename, open="rt")
comments <- readLines(dnds_file, 1) # Read one line 
close(dnds_file)


#' **`r comments`**
#' -----------------------------
#' 
#' 
actual_dnds <- read.table(actual_dnds_filename, header=TRUE, na.strings="None", comment.char = "#", sep=",")
dim(actual_dnds)
head(actual_dnds)
tail(actual_dnds)
str(actual_dnds)
summary(actual_dnds)


expected_dnds <- read.table(expected_dnds_filename, header=TRUE, sep=',')
dim(expected_dnds)
head(expected_dnds)
str(expected_dnds)
summary(expected_dnds)

all(expected_dnds$Site == actual_dnds$Site)
actual_dnds <- actual_dnds[START_CODON:END_CODON,]
summary(actual_dnds)
expected_dnds <- expected_dnds[START_CODON:END_CODON,]
summary(expected_dnds)


# check consistency
all(expected_dnds$Site == actual_dnds$Site)


MIN_SITE <- min(expected_dnds$Site)
MAX_SITE <- max(expected_dnds$Site)
expected_dnds$MultisiteAveDnDs <- apply(expected_dnds, 1, function(row) {
  site <- as.numeric(row["Site"])
  total_dnds <- 0
  total_sites <- 0
  for (adjsite in max(MIN_SITE, site-SMOOTH_DIST):min(MAX_SITE, site+SMOOTH_DIST)) {
    total_dnds <- total_dnds + expected_dnds[expected_dnds$Site==adjsite, "Omega"]
    total_sites <- total_sites +1
  }
  multisite_ave <- total_dnds / total_sites
  return (multisite_ave)
})


#' **Scatterplot actual vs expected dn ds together**

fullDat <- data.frame(actual_dnds, Expected=expected_dnds$Omega, ExpectedMultisite=expected_dnds$MultisiteAveDnDs)
fullDat <- na.omit(fullDat)
summary(fullDat)

ggplot(fullDat[!is.na(fullDat$aveDnDs),], aes(x=aveDnDs, y=Expected)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  xlab("\nInferred") + 
  ylab("Expected\n") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  scale_x_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDat$aveDnDs, fullDat$Expected))) +
  scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDat$aveDnDs, fullDat$Expected))) +
  coord_fixed(ratio=1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Weighted by Substitutions")

ggplot(fullDat[!is.na(fullDat$simpleDnDs),], aes(x=simpleDnDs, y=Expected)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  xlab("\nInferred") + 
  ylab("Expected\n") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  scale_x_continuous(breaks=seq(0, max(fullDat$simpleDnDs), 0.5), limits=c(0, max(fullDat$simpleDnDs, fullDat$Expected))) +
  scale_y_continuous(breaks=seq(0, max(fullDat$simpleDnDs), 0.5), limits=c(0, max(fullDat$simpleDnDs, fullDat$Expected))) +
  coord_fixed(ratio=1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Weighted by Substitutions")

ggplot(fullDat[!is.na(fullDat$dNdSWeightBySubst),], aes(x=dNdSWeightBySubst, y=Expected)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  xlab("\nInferred") + 
  ylab("Expected\n") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  scale_x_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDat$dNdSWeightBySubst, fullDat$Expected))) +
  scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDat$dNdSWeightBySubst, fullDat$Expected))) +
  coord_fixed(ratio=1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Weighted by Substitutions")

ggplot(fullDat[!is.na(fullDat$dNdSWeightByReads),], aes(x=dNdSWeightByReads, y=Expected)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  xlab("\nInferred") + 
  ylab("Expected\n") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  scale_x_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDat$dNdSWeightByReads, fullDat$Expected))) +
  scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDat$dNdSWeightByReads, fullDat$Expected))) +
  coord_fixed(ratio=1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Weighted by Substitutions")

ggplot(fullDat, aes(x=multisiteAvedNdS, y=ExpectedMultisite)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  xlab("\nInferred") + 
  ylab("Expected\n") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  coord_fixed(ratio=1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Averaged Across Surrounding Sites")

ggplot(fullDat, aes(x=multisitedNdSWeightBySubst, y=ExpectedMultisite)) + 
  geom_point() +
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  xlab("\nInferred") + 
  ylab("Expected\n") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  coord_fixed(ratio=1) + 
  ggtitle("Inferred vs Expected dn/ds, Inferred dn/ds Averaged Across Surrounding Sites")


#' **Smoothed Scatterplot of Site dn/ds across the genome**
#+ fig.width=28
fullDatBySource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars=c("Ref", "Site"),
                                              measure.vars=c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", 
                                                             "multisiteAvedNdS", "multisitedNdSWeightBySubst", 
                                                             "Expected", "ExpectedMultisite", "simpleDnDs"),
                                              variable.name="source", value.name="dnds")
head(fullDatBySource)
tail(fullDatBySource)
str(fullDatBySource)
summary(fullDatBySource)
#ggplot(fullDatBySource[fullDatBySource$source %in% c("aveDnDs", "simpleDnDs", "Expected"),], 
ggplot(fullDatBySource[fullDatBySource$source %in% c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "simpleDnDs", "Expected"),], 
       aes(x=Site, y=dnds, color=source, size=2) ) + 
  geom_smooth() + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
  theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
      legend.text=element_text(size=24), legend.title=element_blank())

#' **Line Plot of Site dn/ds across the genome**
#' 
#+ fig.width=20
ggplot(fullDatBySource[fullDatBySource$source %in% c("aveDnDs", "dNdSWeightBySubst", "dNdSWeightByReads", "simpleDnDs", "Expected"),], 
       aes(x=Site, y=dnds, color=source) ) + 
  geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site") + 
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

#' **Plot the expected mutation rate across the genome**
#' 
#+ fig.width=20
ggplot(expected_dnds, aes(x=Site, y=Scaling_factor*BASE_MUTATION_RATE) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Mutation Rate Scaling Factor") + 
  scale_y_continuous(breaks=seq(0, 0.25, 0.01), limits=c(0, 0.25)) +
  ggtitle("Mutation Rate Along Genome")

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

#' **Find how correlated the actual dn-ds vs expected dn/ds are**
dnds_cor <- cor(fullDat$dN_minus_dS, fullDat$Expected, method="spearman", use="complete.obs")
print(dnds_cor)

#' **Find how correlated the actual dn/ds weighted by reads vs expected dn/ds are**
dnds_cor <- cor(fullDat$dNdSWeightByReads, fullDat$Expected, method="spearman", use="complete.obs")
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

#' **Find how concordance correlated the actual dn-ds vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$dN_minus_dS, fullDat$Expected)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual dn/ds weighted by reads vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$dNdSWeightByReads, fullDat$Expected)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual dn/ds smoothed over surrounding sites vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$multisiteAvedNdS, fullDat$ExpectedMultisite)
print(dnds_ccc$rho.c)

#' **Find how concordance correlated the actual dn/ds smoothed over surrounding sites weighted by substitutions vs expected dn/ds are**
dnds_ccc <-  epi.ccc(fullDat$multisitedNdSWeightBySubst, fullDat$ExpectedMultisite)
print(dnds_ccc$rho.c)
