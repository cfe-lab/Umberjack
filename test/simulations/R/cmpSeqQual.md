

```r
library(limma)
library(ggplot2)
library(reshape2)

NUC_PER_CODON <- 3
WINDOWSIZE <- 400
REF_LEN_NUC <- 9000
REF_LEN_CODON <- REF_LEN_NUC%/%NUC_PER_CODON

INDELIBLE_RATES_OUTPUT_COMMENT_LINES <- 10
POPN_SIZE <- 10000
START_NUCPOS <- 3001
END_NUCPOS <- 3400

START_CODON <- 1001  #((START_NUCPOS - 1) / NUC_PER_CODON) + 1
END_CODON <- 1333

actual_dnds_filename <- "../data/out/consensus/mut100_1x_errfree_window400/actual_dnds_by_site.tsv"
dnds_file <- file(actual_dnds_filename, open = "rt")
comments <- readLines(dnds_file, 1)  # Read one line 
close(dnds_file)
```


**# ref=consensus,ref_len=9000,sam=/home/thuy/gitrepo/SlidingWindow/test/simulations/data/indelible/cov1x/sample_genomes.100.1x.errfree.consensus.sam,mapping qual cutoff=19,read qual cutoff=19,max fraction N=0.1,start nuc pos=3001,end nuc pos=3004,windowsize=400,window depth thresh=50,window breadth fraction=0.9,pvalue=0.05**
-----------------------------




```r
actual_dnds <- read.table(actual_dnds_filename, header = TRUE, na.strings = "None", 
    comment.char = "#")
dim(actual_dnds)
```

```
## [1] 3000    9
```

```r
head(actual_dnds)
```

```
##         Ref Site dNdS dN_minus_dS Windows Codons NonSyn Syn Subst
## 1 consensus    1   NA          NA       0     NA     NA  NA    NA
## 2 consensus    2   NA          NA       0     NA     NA  NA    NA
## 3 consensus    3   NA          NA       0     NA     NA  NA    NA
## 4 consensus    4   NA          NA       0     NA     NA  NA    NA
## 5 consensus    5   NA          NA       0     NA     NA  NA    NA
## 6 consensus    6   NA          NA       0     NA     NA  NA    NA
```

```r
tail(actual_dnds)
```

```
##            Ref Site dNdS dN_minus_dS Windows Codons NonSyn Syn Subst
## 2995 consensus 2995   NA          NA       0     NA     NA  NA    NA
## 2996 consensus 2996   NA          NA       0     NA     NA  NA    NA
## 2997 consensus 2997   NA          NA       0     NA     NA  NA    NA
## 2998 consensus 2998   NA          NA       0     NA     NA  NA    NA
## 2999 consensus 2999   NA          NA       0     NA     NA  NA    NA
## 3000 consensus 3000   NA          NA       0     NA     NA  NA    NA
```

```r
str(actual_dnds)
```

```
## 'data.frame':	3000 obs. of  9 variables:
##  $ Ref        : Factor w/ 1 level "consensus": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Site       : int  1 2 3 4 5 6 7 8 9 10 ...
##  $ dNdS       : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ dN_minus_dS: num  NA NA NA NA NA NA NA NA NA NA ...
##  $ Windows    : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ Codons     : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ NonSyn     : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ Syn        : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ Subst      : num  NA NA NA NA NA NA NA NA NA NA ...
```

```r
summary(actual_dnds)
```

```
##         Ref            Site           dNdS       dN_minus_dS  
##  consensus:3000   Min.   :   1   Min.   :0.0    Min.   :-3.4  
##                   1st Qu.: 751   1st Qu.:0.1    1st Qu.:-1.4  
##                   Median :1500   Median :0.3    Median :-1.1  
##                   Mean   :1500   Mean   :0.4    Mean   :-1.1  
##                   3rd Qu.:2250   3rd Qu.:0.5    3rd Qu.:-0.9  
##                   Max.   :3000   Max.   :2.4    Max.   : 2.1  
##                                  NA's   :2901   NA's   :2901  
##     Windows           Codons         NonSyn           Syn      
##  Min.   :0.0000   Min.   :1795   Min.   :  1.0   Min.   : 8.5  
##  1st Qu.:0.0000   1st Qu.:3072   1st Qu.: 11.1   1st Qu.:30.1  
##  Median :0.0000   Median :3078   Median : 25.0   Median :34.0  
##  Mean   :0.0657   Mean   :2949   Mean   : 30.5   Mean   :37.2  
##  3rd Qu.:0.0000   3rd Qu.:3092   3rd Qu.: 36.9   3rd Qu.:40.4  
##  Max.   :2.0000   Max.   :3092   Max.   :158.6   Max.   :83.4  
##                   NA's   :2901   NA's   :2901    NA's   :2901  
##      Subst      
##  Min.   : 27.4  
##  1st Qu.: 48.0  
##  Median : 61.1  
##  Mean   : 67.7  
##  3rd Qu.: 74.9  
##  Max.   :184.1  
##  NA's   :2901
```

```r


expected_dnds <- read.table("../data/indelible/sample_genomes.100.rates.csv", 
    header = TRUE, sep = ",")
dim(expected_dnds)
```

```
## [1] 3000    5
```

```r
head(expected_dnds)
```

```
##   Site Interval Scaling_factor Rate_class Omega
## 1    1        0            100          5  0.55
## 2    2        0            100          1  0.15
## 3    3        0            100         20  2.05
## 4    4        0            100         14  1.45
## 5    5        0            100          2  0.25
## 6    6        0            100          0  0.05
```

```r
str(expected_dnds)
```

```
## 'data.frame':	3000 obs. of  5 variables:
##  $ Site          : int  1 2 3 4 5 6 7 8 9 10 ...
##  $ Interval      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ Scaling_factor: num  100 100 100 100 100 100 100 100 100 100 ...
##  $ Rate_class    : int  5 1 20 14 2 0 1 2 3 3 ...
##  $ Omega         : num  0.55 0.15 2.05 1.45 0.25 0.05 0.15 0.25 0.35 0.35 ...
```

```r
summary(expected_dnds)
```

```
##       Site         Interval Scaling_factor   Rate_class        Omega      
##  Min.   :   1   Min.   :0   Min.   :100    Min.   : 0.00   Min.   :0.050  
##  1st Qu.: 751   1st Qu.:0   1st Qu.:100    1st Qu.: 1.00   1st Qu.:0.150  
##  Median :1500   Median :0   Median :100    Median : 3.00   Median :0.350  
##  Mean   :1500   Mean   :0   Mean   :100    Mean   : 4.56   Mean   :0.506  
##  3rd Qu.:2250   3rd Qu.:0   3rd Qu.:100    3rd Qu.: 6.00   3rd Qu.:0.650  
##  Max.   :3000   Max.   :0   Max.   :100    Max.   :33.00   Max.   :3.350
```

```r

actual_dnds <- actual_dnds[1001:1133, ]
summary(actual_dnds)
```

```
##         Ref           Site           dNdS       dN_minus_dS   
##  consensus:133   Min.   :1001   Min.   :0.01   Min.   :-3.43  
##                  1st Qu.:1034   1st Qu.:0.14   1st Qu.:-1.39  
##                  Median :1067   Median :0.27   Median :-1.15  
##                  Mean   :1067   Mean   :0.41   Mean   :-1.05  
##                  3rd Qu.:1100   3rd Qu.:0.48   3rd Qu.:-0.86  
##                  Max.   :1133   Max.   :2.38   Max.   : 2.05  
##                                 NA's   :34     NA's   :34     
##     Windows         Codons         NonSyn           Syn      
##  Min.   :0.00   Min.   :1795   Min.   :  1.0   Min.   : 8.5  
##  1st Qu.:0.00   1st Qu.:3072   1st Qu.: 11.1   1st Qu.:30.1  
##  Median :2.00   Median :3078   Median : 25.0   Median :34.0  
##  Mean   :1.48   Mean   :2949   Mean   : 30.5   Mean   :37.2  
##  3rd Qu.:2.00   3rd Qu.:3092   3rd Qu.: 36.9   3rd Qu.:40.4  
##  Max.   :2.00   Max.   :3092   Max.   :158.6   Max.   :83.4  
##                 NA's   :34     NA's   :34      NA's   :34    
##      Subst      
##  Min.   : 27.4  
##  1st Qu.: 48.0  
##  Median : 61.1  
##  Mean   : 67.7  
##  3rd Qu.: 74.9  
##  Max.   :184.1  
##  NA's   :34
```

```r
expected_dnds <- expected_dnds[1001:1133, ]
summary(expected_dnds)
```

```
##       Site         Interval Scaling_factor   Rate_class       Omega      
##  Min.   :1001   Min.   :0   Min.   :100    Min.   : 0.0   Min.   :0.050  
##  1st Qu.:1034   1st Qu.:0   1st Qu.:100    1st Qu.: 2.0   1st Qu.:0.250  
##  Median :1067   Median :0   Median :100    Median : 3.0   Median :0.350  
##  Mean   :1067   Mean   :0   Mean   :100    Mean   : 4.2   Mean   :0.469  
##  3rd Qu.:1100   3rd Qu.:0   3rd Qu.:100    3rd Qu.: 6.0   3rd Qu.:0.650  
##  Max.   :1133   Max.   :0   Max.   :100    Max.   :17.0   Max.   :1.750
```


**Scatterplot actual vs expected dn ds together**


```r
# check consistency
all(expected_dnds$Site == actual_dnds$Site)
```

```
## [1] TRUE
```

```r

fullDat <- data.frame(expected_site = expected_dnds$Site, expected = expected_dnds$Omega, 
    actual_site = actual_dnds$Site, actual = actual_dnds$dNdS)
head(fullDat[!is.na(fullDat$actual), ])
```

```
##    expected_site expected actual_site actual
## 1           1001     0.15        1001 0.3808
## 2           1002     0.15        1002 0.2031
## 3           1003     1.45        1003 2.3788
## 4           1004     0.25        1004 0.3695
## 8           1008     1.35        1008 1.5057
## 10          1010     0.15        1010 0.2708
```

```r
ggplot(fullDat, aes(x = actual, y = expected)) + geom_smooth(method = lm)
```

```
## Warning: Removed 34 rows containing missing values (stat_smooth).
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


**Scatterplot the dn/ds across the genome**


```r
fullDatBySource <- reshape2:::melt.data.frame(data = fullDat, na.rm = FALSE, 
    id.vars = "actual_site", measure.vars = c("expected", "actual"), variable.name = "source", 
    value.name = "dnds")
head(fullDatBySource)
```

```
##   actual_site   source dnds
## 1        1001 expected 0.15
## 2        1002 expected 0.15
## 3        1003 expected 1.45
## 4        1004 expected 0.25
## 5        1005 expected 0.75
## 6        1006 expected 0.35
```

```r
tail(fullDatBySource)
```

```
##     actual_site source   dnds
## 261        1128 actual 0.4519
## 262        1129 actual 0.1843
## 263        1130 actual 0.2549
## 264        1131 actual 0.5977
## 265        1132 actual 0.1027
## 266        1133 actual 0.2668
```

```r
str(fullDatBySource)
```

```
## 'data.frame':	266 obs. of  3 variables:
##  $ actual_site: int  1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 ...
##  $ source     : Factor w/ 2 levels "expected","actual": 1 1 1 1 1 1 1 1 1 1 ...
##  $ dnds       : num  0.15 0.15 1.45 0.25 0.75 0.35 1.05 1.35 0.55 0.15 ...
```

```r
summary(fullDatBySource)
```

```
##   actual_site        source         dnds     
##  Min.   :1001   expected:133   Min.   :0.01  
##  1st Qu.:1034   actual  :133   1st Qu.:0.15  
##  Median :1067                  Median :0.35  
##  Mean   :1067                  Mean   :0.44  
##  3rd Qu.:1100                  3rd Qu.:0.55  
##  Max.   :1133                  Max.   :2.38  
##                                NA's   :34
```

```r
ggplot(fullDatBySource, aes(x = site, y = dnds, color = source)) + geom_smooth() + 
    xlab("Codon Site Along Genome") + ylab("dN/dS") + ggtitle("dn/ds by site")
```

```
## Error: object 'site' not found
```

```r

ggplot(fullDatBySource, aes(x = actual_site, y = dnds, color = source)) + geom_line() + 
    xlab("Codon Site Along Genome") + ylab("dN/dS") + ggtitle("dn/ds by site")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r


ggplot(fullDatBySource, aes(x = actual_site, y = dnds, color = source)) + geom_point() + 
    xlab("Codon Site Along Genome") + ylab("dN/dS") + ggtitle("dn/ds by site")
```

```
## Warning: Removed 34 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 


**Plot the unambiguous codon depth across genome**


```r
ggplot(actual_dnds, aes(x = Site, y = Codons)) + geom_line() + xlab("Codon Site Along Genome") + 
    ylab("Sequences with Unambiguous Codons") + ggtitle("Population Unambiguous Codons Across Genome")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


**Plot the nonsynonymous substitutions across genome**


```r
ggplot(actual_dnds, aes(x = Site, y = NonSyn)) + geom_line() + xlab("Codon Site Along Genome") + 
    ylab("Nonsynonymous Substitutions") + ggtitle("Population Nonsynonymous Substitutions Across Genome")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


**Plot the synonymous substitutions across genome**


```r
ggplot(actual_dnds, aes(x = Site, y = Syn)) + geom_line() + xlab("Codon Site Along Genome") + 
    ylab("Synonymous Substitutions") + ggtitle("Population Synonymous Substitutions Across Genome")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


**Plot the substitutions across genome**


```r
ggplot(actual_dnds, aes(x = Site, y = Subst)) + geom_line() + xlab("Codon Site Along Genome") + 
    ylab("Substitutions") + ggtitle("Population Substitutions Across Genome")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


**Plot the Windows across genome**


```r
ggplot(actual_dnds, aes(x = Site, y = Windows)) + geom_line() + xlab("Codon Site Along Genome") + 
    ylab("Windows") + ggtitle("Windows Across Genome")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


**Plot the expected mutation rate across the genome**


```r
ggplot(expected_dnds, aes(x = Site, y = Scaling_factor)) + geom_line() + xlab("Codon Site Along Genome") + 
    ylab("Mutation Rate Scaling Factor") + ggtitle("Mutation Along Genome")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


**Plot the Expected Omega rate across the genome**


```r
ggplot(expected_dnds, aes(x = Site, y = Omega)) + geom_line() + xlab("Codon Site Along Genome") + 
    ylab("dn/dS Expected") + ggtitle("Expected Selection Along Genome")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 


**Paired test without assumption of normalcy**


```r
htest <- wilcox.test(actual_dnds$dNdS, expected_dnds$Omega, paired = TRUE, exact = TRUE, 
    na.action = "na.exclude")
print(htest)
```

```
## 
## 	Wilcoxon signed rank test
## 
## data:  actual_dnds$dNdS and expected_dnds$Omega
## V = 2681, p-value = 0.4749
## alternative hypothesis: true location shift is not equal to 0
```

```r
# dnds_cor <- cor(log(actual_dnds$dNdS), expected_dnds$Omega,
# method='spearman', use='pairwise.complete.obs')
dnds_cor <- cor(actual_dnds$dNdS, expected_dnds$Omega, method = "spearman", 
    use = "complete.obs")
print(dnds_cor)
```

```
## [1] 0.9093
```

```r
dnds_dnMinusds_cor <- cor(actual_dnds$dN_minus_dS, expected_dnds$Omega, method = "spearman", 
    use = "complete.obs")
print(dnds_dnMinusds_cor)
```

```
## [1] 0.5903
```


**Linear Regression to see how dN/dS varies by reads, substitutions, seq error**



```r
# rows are codon sites. cols are samples our samples are actual, expected
# linDat <- data.frame(row.names=expected_dnds$Site,
# site=expected_dnds$Site, expected_dnds=expected_dnds$Omega,
# actual_dnds=actual_dnds$dNdS ) head(linDat[!is.na(linDat$actual),])
# str(linDat) summary(linDat) linDes <-
# reshape2:::melt.data.frame(data=linDat, na.rm = FALSE, id.vars='site',
# measure.vars=c('expected_dnds', 'actual_dnds'), variable.name='source',
# value.name='dnds') linDesModMat <- model.matrix(~source, linDes)
# str(linDesModMat) fit <- lmFit(linDat, linDesModMat) ebfit <- eBayes(fit)
# tophits <- topTable(wtEbFit, coef = grep('devStage',
# colnames(coef(wtEbFit)))))
```


