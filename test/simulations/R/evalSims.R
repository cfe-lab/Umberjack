library(limma)
library(ggplot2)
library(reshape2)

INDELIBLE_RATES_OUTPUT_COMMENT_LINES <- 10
ACTUAL_DNDS_COMMENT_LINES <- 1



dnds_file <- file("../data/out/consensus/actual_dnds_by_site.tsv", open="rt")
comments <- readLines(dnds_file, 1) # Read one line 

#' **`r comments`**
#' -----------------------------
#' 
#' 
actual_dnds <- read.table("../data/out/consensus/actual_dnds_by_site.tsv", header=TRUE, na.strings="None", comment.char = "#")
dim(actual_dnds)
head(actual_dnds)
tail(actual_dnds)
str(actual_dnds)
summary(actual_dnds)
# Convert NA Codons to zero
actual_dnds$Codons[is.na(actual_dnds$Codons)] <- 0
actual_dnds$NonSyn[is.na(actual_dnds$NonSyn)] <- 0
actual_dnds$Syn[is.na(actual_dnds$Syn)] <- 0
actual_dnds$Subst[is.na(actual_dnds$Subst)] <- 0
summary(actual_dnds)

expected_dnds <- read.table("../data/sample_genomes.rates", header=TRUE, sep=',')
dim(expected_dnds)
head(expected_dnds)
str(expected_dnds)
summary(expected_dnds)


# orig_dnds <- read.table("../data/sample_genomes.rates.orig", header=TRUE, sep=',')
# dim(orig_dnds)
# head(orig_dnds)
# str(orig_dnds)
# summary(orig_dnds)
# all(orig_dnds$omega == expected_dnds$Omega)
# sum(orig_dnds$omega == expected_dnds$Omega)

#' **Paired test without assumption of normalcy**
htest <-  wilcox.test(actual_dnds$dNdS, expected_dnds$Omega, paired=TRUE, exact=TRUE, na.action="na.exclude")
print (htest)

#' **Scatterplot actual vs expected dn ds together**
#+ fig.width=24
fullDat <- data.frame(expected_site=expected_dnds$Site,
                      expected=expected_dnds$Omega,
                      actual_site=actual_dnds$Site,
                      actual=actual_dnds$dNdS)
head(fullDat[!is.na(fullDat$actual),])
ggplot(fullDat, aes(x=actual, y=expected)) + geom_smooth(method=lm)

#' **Scatterplot the dn/ds across the genome**
#+ fig.width=28
fullDatBySource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars="actual_site", measure.vars=c("expected", "actual"),
                                              variable.name="source", value.name="dnds")
head(fullDatBySource)
tail(fullDatBySource)
str(fullDatBySource)
summary(fullDatBySource)
ggplot(fullDatBySource, aes(x=site, y=dnds, color=source) ) + geom_smooth() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site")

ggplot(fullDatBySource, aes(x=actual_site, y=dnds, color=source) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site")


ggplot(fullDatBySource, aes(x=actual_site, y=dnds, color=source) ) + geom_point() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site")


#' **Plot the unambiguous codon depth across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=Codons) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Sequences with Unambiguous Codons") + 
  ggtitle("Population Unambiguous Codons Across Genome")


#' **Plot the nonsynonymous substitutions across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=NonSyn) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Nonsynonymous Substitutions") + 
  ggtitle("Population Nonsynonymous Substitutions Across Genome")


#' **Plot the synonymous substitutions across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=Syn) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Synonymous Substitutions") + 
  ggtitle("Population Synonymous Substitutions Across Genome")


#' **Plot the substitutions across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=Subst) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Substitutions") + 
  ggtitle("Population Substitutions Across Genome")

#' **Plot the Windows across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=Windows) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Windows") + 
  ggtitle("Windows Across Genome")

#' **Plot the expected mutation rate across the genome**
#+ fig.width=24
ggplot(expected_dnds, aes(x=Site, y=Scaling_factor) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Mutation Rate Scaling Factor") + 
  ggtitle("Mutation Along Genome")



#' **Plot the Expected Omega rate across the genome**
#+ fig.width=24
ggplot(expected_dnds, aes(x=Site, y=Omega) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dn/dS Expected") + 
  ggtitle("Expected Selection Along Genome")

#+ **Find how correlated the actual vs expected dn/ds are**
#dnds_cor <- cor(log(actual_dnds$dNdS), expected_dnds$Omega, method="spearman", use="pairwise.complete.obs")
dnds_cor <- cor(actual_dnds$dNdS, expected_dnds$Omega, method="spearman", use="complete.obs")
print(dnds_cor)

#+ **Find how correlated the actual dn-ds vs expected dn/ds are**
dnds_dnMinusds_cor <- cor(actual_dnds$dN_minus_dS, expected_dnds$Omega, method="spearman", use="complete.obs")
print(dnds_dnMinusds_cor)

