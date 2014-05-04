# Double check data integrity

expected_dnds <- read.table("../data/sample_genomes.rates", header=TRUE, sep=',')
dim(expected_dnds)
head(expected_dnds)
str(expected_dnds)
summary(expected_dnds)

# Created by Art originally
orig_dnds <- read.table("../data/sample_genomes.rates.orig", header=TRUE, sep=',')
dim(orig_dnds)
head(orig_dnds)
str(orig_dnds)
summary(orig_dnds)
all(orig_dnds$omega == expected_dnds$Omega)
sum(orig_dnds$omega == expected_dnds$Omega)

# Check that INDElible uses the same codon substitution rate classes in the same site for every tree it generated
INDELIBLE_RATE_COMMENT_LINES <- 9
scaling_factors <- c("1.0", "2.0", "5.0", "10.0", "20.0", "50.0", "100.0")


rates_filename <- paste("../data/scaling_",  scaling_factors[1], "_RATES.txt", sep="")
print (rates_filename)
indelible_rates1 <- read.table(rates_filename, header=TRUE, skip=INDELIBLE_RATE_COMMENT_LINES, fill=TRUE)
dim(indelible_rates1)
head(indelible_rates1)
str(indelible_rates1)
summary(indelible_rates1)

checkIndelibleRates <- function (scaling_factor) {
  rates_filename <- paste("../data/scaling_",  scaling_factor, "_RATES.txt", sep="")
  print (rates_filename)
  indelible_rates <- read.table(rates_filename, header=TRUE, skip=INDELIBLE_RATE_COMMENT_LINES, fill=TRUE)
  print(dim(indelible_rates))
  print(head(indelible_rates))
  print(str(indelible_rates))
  print(summary(indelible_rates))
  result <- all(indelible_rates1 == indelible_rates)
  return (result)
}

indelibleRatesSameAs1 <- lapply(scaling_factors, checkIndelibleRates)
print(indelibleRatesSameAs1)



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



actual_dnds <- read.table("../data/scaling_100.0_TRUE.3001_3400.dnds.nobranchcorr.R.tsv", header=TRUE, na.strings="None", comment.char = "#")
dim(actual_dnds)
head(actual_dnds)
tail(actual_dnds)
str(actual_dnds)
summary(actual_dnds)
# convert dn/ds with bad pvalues to NA
actual_dnds$NegPvalue <- as.numeric(as.character(actual_dnds$NegPvalue))
actual_dnds$dNdS[actual_dnds$Pvalue>0.05 & actual_dnds$NegPvalue>0.05] <- NA
actual_dnds$ScaledDnMinusDs[actual_dnds$Pvalue>0.05 & actual_dnds$NegPvalue>0.05] <- NA
summary(actual_dnds)
# Convert NA to zero
# actual_dnds$dNdS[is.na(actual_dnds$dNdS)] <- 0
# actual_dnds$ScaledDnMinusDs[is.na(actual_dnds$ScaledDnMinusDs)] <- 0
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
htest <-  wilcox.test(actual_dnds$dNdS, expected_dnds$Omega[1001:1133], paired=TRUE, exact=TRUE, na.action="na.exclude")
print (htest)

#' **Scatterplot actual vs expected dn ds together**
#+ fig.width=24
fullDat <- data.frame(expected_site=expected_dnds$Site[1001:1133],
                      expected=expected_dnds$Omega[1001:1133],
                      actual_site=actual_dnds$Site,
                      actual=actual_dnds$dNdS)
                      
                      
head(fullDat[!is.na(fullDat$actual),])
ggplot(fullDat, aes(x=actual, y=expected)) + geom_smooth(method=lm)

#' **Scatterplot the dn/ds across the genome**
#+ fig.width=28
fullDatBySource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars="actual_site", measure.vars=c("actual", "expected"), variable.name="source", value.name="dnds")
head(fullDatBySource)
tail(fullDatBySource)
str(fullDatBySource)
summary(fullDatBySource)
ggplot(fullDatBySource, aes(x=actual_site, y=dnds, color=source) ) + geom_smooth() + 
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






#' **Plot the nonsynonymous substitutions across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=NS) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Nonsynonymous Substitutions") + 
  ggtitle("Population Nonsynonymous Substitutions Across Genome")


#' **Plot the synonymous substitutions across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=S) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Synonymous Substitutions") + 
  ggtitle("Population Synonymous Substitutions Across Genome")


#' **Plot the substitutions across genome**
#+ fig.width=24
ggplot(actual_dnds, aes(x=Site, y=NS + S) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("Substitutions") + 
  ggtitle("Population Substitutions Across Genome")


#' **Plot the Expected Omega rate across the genome**
#+ fig.width=24
# ggplot(expected_dnds, aes(x=Site, y=Omega[927:1059]) ) + geom_line() + 
#   xlab("Codon Site Along Genome") + 
#   ylab("dn/dS Expected") + 
#   ggtitle("Expected Selection Along Genome")
ggplot(expected_dnds[1001:1133,], aes(x=Site, y=Omega) ) + geom_line() + 
  xlab("Codon Site Along Genome") + 
  ylab("dn/dS Expected") + 
  ggtitle("Expected Selection Along Genome")


 
#+ **Find how correlated the actual vs expected dn/ds are**
#dnds_cor <- cor(log(actual_dnds$dNdS), expected_dnds$Omega, method="spearman", use="pairwise.complete.obs")
dnds_cor <- cor(actual_dnds$dNdS, expected_dnds$Omega[1001:1133], method="spearman", use="complete.obs")
print(dnds_cor)


#+ **Find how correlated the actual dn-ds vs expected dn/ds are**
dnds_dnMinusds_cor <- cor(actual_dnds$dN_minus_dS, expected_dnds$Omega[1001:1133], method="spearman", use="complete.obs")
print(dnds_dnMinusds_cor)

