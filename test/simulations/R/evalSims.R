library(limma)
library(ggplot2)

actual_dnds <- read.table("../data/actual_dnds.tsv", header=TRUE)
expected_dnds <- read.table("../data/expected_dnds.tsv", header=TRUE)  # TODO:  make this file

# paired test without assumption of normalcy
htest <-  wilcox.test(actual_dnds, expected_dnds, paired=TRUE, exact=TRUE, conf.int=TRUE, conf.level=0.95)
print (htest$statistic)
print (htest$p.value)
print (htest$conf.int)


# scatterplot actual vs expected dn ds together
fullDat <- data.frame(actual=actual_dnds$dnds, expected=expected_dnds$dnds)
ggplot(fullDat, aes(x=actual, y=expected)) + geom_smooth(method=lm)

# Find how correlated the actual vs expected dnds are
dnds_cor <- cor(actual_dnds, expected_dnds)
heatmap(dnds_cor)