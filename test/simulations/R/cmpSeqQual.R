library(limma)
library(ggplot2)
library(reshape2)

NUC_PER_CODON <- 3
WINDOWSIZE <- 400
REF_LEN_NUC <- 9000
REF_LEN_CODON <- REF_LEN_NUC %/% NUC_PER_CODON

INDELIBLE_RATES_OUTPUT_COMMENT_LINES <- 10
POPN_SIZE <- 10000
START_NUCPOS <- 3001
END_NUCPOS <- 3400

START_CODON <- 1001 #((START_NUCPOS - 1) / NUC_PER_CODON) + 1
END_CODON <- 1333
  
#actual_dnds_filename <- "../data/out/consensus/mut100_1x_errfree_window400/actual_dnds_by_site.tsv"
actual_dnds_filename <- "../data/out/consensus/mut100_1x_window400/actual_dnds_by_site.tsv"
dnds_file <- file(actual_dnds_filename, open="rt")
comments <- readLines(dnds_file, 1) # Read one line 
close(dnds_file)

#' **`r comments`**
#' -----------------------------
#' 
#' 
actual_dnds <- read.table(actual_dnds_filename, header=TRUE, na.strings="None", comment.char = "#")
dim(actual_dnds)
head(actual_dnds)
tail(actual_dnds)
str(actual_dnds)
summary(actual_dnds)


expected_dnds <- read.table("../data/indelible/sample_genomes.100.rates.csv", header=TRUE, sep=',')
dim(expected_dnds)
head(expected_dnds)
str(expected_dnds)
summary(expected_dnds)

all(expected_dnds$Site == actual_dnds$Site)
actual_dnds <- actual_dnds[1001:1133,]
summary(actual_dnds)
expected_dnds <- expected_dnds[1001:1133,]
summary(expected_dnds)
#' **Scatterplot actual vs expected dn ds together**
#+ fig.width=24
# check consistency
all(expected_dnds$Site == actual_dnds$Site)

fullDat <- data.frame(Site=expected_dnds$Site,
                      Expected=expected_dnds$Omega,
                      Inferred=actual_dnds$dNdS)
head(fullDat[!is.na(fullDat$Inferred),])
fullDatNoNA <- fullDat[!is.na(fullDat$Inferred),]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scale_fill_manual(values=cbPalette)
#ggplot(fullDatNoNA, aes(x=Inferred, y=Expected)) + geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
ggplot(fullDatNoNA, aes(x=Inferred, y=Expected)) + geom_smooth(method=lm, size=4, color="#006600", fill="#00FF99") +
  xlab("\nInferred") + 
  ylab("Expected\n") + 
  theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
  scale_x_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDatNoNA$Inferred))) +
  scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, max(fullDatNoNA$Expected))) +
  coord_fixed(ratio=1)

#' **Scatterplot the dn/ds across the genome**
#+ fig.width=28
fullDatBySource <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars="site", measure.vars=c("Expected", "Inferred"),
                                              variable.name="source", value.name="dnds")
head(fullDatBySource)
tail(fullDatBySource)
str(fullDatBySource)
summary(fullDatBySource)
ggplot(fullDatBySource, aes(x=site, y=dnds, color=source) ) + geom_smooth() + 
  xlab("Codon Site") + 
  ylab("dN/dS") + 
  #ggtitle("dn/ds by site")
theme(plot.title=element_text(size=36), axis.title=element_text(size=32), axis.text=element_text(size=24), 
      legend.text=element_text(size=24), legend.title=element_blank())

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

#' **Paired test without assumption of normalcy**
htest <-  wilcox.test(actual_dnds$dNdS, expected_dnds$Omega, paired=TRUE, exact=TRUE, na.action="na.exclude")
print (htest)

#+ **Find how correlated the actual vs expected dn/ds are**
#dnds_cor <- cor(log(actual_dnds$dNdS), expected_dnds$Omega, method="spearman", use="pairwise.complete.obs")
dnds_cor <- cor(actual_dnds$dNdS, expected_dnds$Omega, method="spearman", use="complete.obs")
print(dnds_cor)

#+ **Find how correlated the actual dn-ds vs expected dn/ds are**
dnds_dnMinusds_cor <- cor(actual_dnds$dN_minus_dS, expected_dnds$Omega, method="spearman", use="complete.obs")
print(dnds_dnMinusds_cor)


#' **Linear Regression to see how dN/dS varies by reads, substitutions, seq error**
#' 

# rows are codon sites. cols are samples  
# our samples are actual, expected  
# linDat <- data.frame(row.names=expected_dnds$Site,
#                      site=expected_dnds$Site,
#                      expected_dnds=expected_dnds$Omega,                     
#                      actual_dnds=actual_dnds$dNdS
#                      )
# head(linDat[!is.na(linDat$actual),])
# str(linDat)
# summary(linDat)
# 
# 
# linDes <- reshape2:::melt.data.frame(data=linDat, na.rm = FALSE, id.vars="site", measure.vars=c("expected_dnds", "actual_dnds"),
#                                      variable.name="source", value.name="dnds")
# linDesModMat <- model.matrix(~source, linDes)
# str(linDesModMat)
# 
# fit <- lmFit(linDat, linDesModMat)
# ebfit <- eBayes(fit)
# 
# tophits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit)))))