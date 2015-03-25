library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)
dnds_csv <- args[1]
dnds_img <- args[2]
# Cols:  Ref  Site  aveDnDs	dNdSWeightBySubst	dN_minus_dS	Windows	Codons	NonSyn	Syn	Subst	dNdSWeightByReads	multisiteAvedNdS	multisitedNdSWeightBySubst	dNdSWeightByReadsNoLowSyn
dndsdat <- read.table(dnds_csv, header=TRUE, comment.char = "#", sep=",")

if (sum(!is.na(dndsdat$dNdSWeightByReadsNoLowSyn)) == 0) {
  warning(paste0("There is no data to plot in ", dnds_csv))
} else {
  fig <- ggplot(dndsdat, aes(x=Site, y=dNdSWeightByReadsNoLowSyn) ) + 
    geom_line(na.rm=TRUE) + 
    xlab("Codon Site") + 
    ylab("dN/dS") + 
    ggtitle("Site Average dn/ds\nExcluding Sites < 1 Substitution")

  ggsave(filename=dnds_img, plot=fig, device=pdf)

}