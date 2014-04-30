library(limma)
library(ggplot2)
library(reshape2)

INDELIBLE_RATES_OUTPUT_COMMENT_LINES <- 10
ACTUAL_DNDS_COMMENT_LINES <- 0

actual_dnds <- read.table("../data/out/consensus/actual_dnds_by_site.tsv", header=TRUE, na.strings="None", skip=ACTUAL_DNDS_COMMENT_LINES)
dim(actual_dnds)
head(actual_dnds)
tail(actual_dnds)
str(actual_dnds)
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
fullDat <- data.frame(site=expected_dnds$Site,
                      actual=actual_dnds$dNdS, 
                      expected=expected_dnds$Omega)
head(fullDat[!is.na(fullDat$actual),])
ggplot(fullDat, aes(x=actual, y=expected)) + geom_smooth(method=lm)

#' **Scatterplot the dn-ds across the genome**
#+ fig.width=16
fullDatBySource <- reshape2:::melt.data.frame(data=fullDat, na.rm = TRUE, id.vars="site", variable.name="source", value.name="dnds")
head(fullDatBySource)
tail(fullDatBySource)
str(fullDatBySource)
summary(fullDatBySource)
ggplot(fullDatBySource, aes(x=site, y=dnds, color=source) ) + geom_point() + 
  xlab("Codon Site Along Genome") + 
  ylab("dN/dS") + 
  ggtitle("dn/ds by site")

# 
#' **Plot the mutation rate across the genome**
#+ fig.width=16
ggplot(expected_dnds, aes(x=Site, y=Scaling_factor) ) + geom_point() + 
  xlab("Codon Site Along Genome") + 
  ylab("Mutation Rate Scaling Factor") + 
  ggtitle("Mutation Along Genome")

# 
#' **Plot the Expected Omega rate across the genome**
#+ fig.width=16
ggplot(expected_dnds, aes(x=Site, y=Omega) ) + geom_point() + 
  xlab("Codon Site Along Genome") + 
  ylab("dn/dS Expected") + 
  ggtitle("Expected Selection Along Genome")

#+ **Find how correlated the actual vs expected dnds are**
#dnds_cor <- cor(log(actual_dnds$dNdS), expected_dnds$Omega, method="spearman", use="pairwise.complete.obs")
dnds_cor <- cor(actual_dnds$dNdS, expected_dnds$Omega, method="spearman", use="complete.obs")
print(dnds_cor)
