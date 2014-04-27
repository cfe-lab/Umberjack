

```r
library(limma)
library(ggplot2)
library(reshape2)

INDELIBLE_RATES_OUTPUT_COMMENT_LINES <- 10

actual_dnds <- read.table("../data/out/consensus/ave_dnds_by_site.tsv", header = TRUE, 
    na.strings = "None")
dim(actual_dnds)
```

```
## [1] 3000    3
```

```r
head(actual_dnds)
```

```
##         Ref Site     dNdS
## 1 consensus    1 4.213074
## 2 consensus    2 0.008495
## 3 consensus    3 0.013434
## 4 consensus    4 0.113133
## 5 consensus    5 0.206868
## 6 consensus    6 0.128783
```

```r
tail(actual_dnds)
```

```
##            Ref Site dNdS
## 2995 consensus 2995   NA
## 2996 consensus 2996   NA
## 2997 consensus 2997   NA
## 2998 consensus 2998   NA
## 2999 consensus 2999   NA
## 3000 consensus 3000   NA
```

```r
str(actual_dnds)
```

```
## 'data.frame':	3000 obs. of  3 variables:
##  $ Ref : Factor w/ 1 level "consensus": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Site: int  1 2 3 4 5 6 7 8 9 10 ...
##  $ dNdS: num  4.21307 0.00849 0.01343 0.11313 0.20687 ...
```

```r
summary(actual_dnds)
```

```
##         Ref            Site           dNdS     
##  consensus:3000   Min.   :   1   Min.   :-0.1  
##                   1st Qu.: 751   1st Qu.: 0.0  
##                   Median :1500   Median : 0.0  
##                   Mean   :1500   Mean   : 0.0  
##                   3rd Qu.:2250   3rd Qu.: 0.0  
##                   Max.   :3000   Max.   : 4.2  
##                                  NA's   :2446
```

```r

expected_dnds <- read.table("../data/sample_genomes.rates", header = TRUE, sep = ",")
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
## 1    1        0              1          3  0.35
## 2    2        0              1          4  0.45
## 3    3        0              1          4  0.45
## 4    4        0              1          7  0.75
## 5    5        0              1          0  0.05
## 6    6        0              1          3  0.35
```

```r
str(expected_dnds)
```

```
## 'data.frame':	3000 obs. of  5 variables:
##  $ Site          : int  1 2 3 4 5 6 7 8 9 10 ...
##  $ Interval      : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ Scaling_factor: num  1 1 1 1 1 1 1 1 1 1 ...
##  $ Rate_class    : int  3 4 4 7 0 3 7 2 3 1 ...
##  $ Omega         : num  0.35 0.45 0.45 0.75 0.05 0.35 0.75 0.25 0.35 0.15 ...
```

```r
summary(expected_dnds)
```

```
##       Site         Interval Scaling_factor    Rate_class   
##  Min.   :   1   Min.   :0   Min.   :  1.0   Min.   : 0.00  
##  1st Qu.: 751   1st Qu.:3   1st Qu.:  5.0   1st Qu.: 1.00  
##  Median :1500   Median :4   Median : 10.0   Median : 3.00  
##  Mean   :1500   Mean   :4   Mean   : 37.7   Mean   : 4.44  
##  3rd Qu.:2250   3rd Qu.:5   3rd Qu.:100.0   3rd Qu.: 6.00  
##  Max.   :3000   Max.   :8   Max.   :100.0   Max.   :27.00  
##      Omega      
##  Min.   :0.050  
##  1st Qu.:0.150  
##  Median :0.350  
##  Mean   :0.494  
##  3rd Qu.:0.650  
##  Max.   :2.750
```


**Paired test without assumption of normalcy**


```r
htest <- wilcox.test(actual_dnds$dNdS, expected_dnds$Omega, paired = TRUE, exact = TRUE, 
    conf.int = TRUE, conf.level = 0.95, na.action = "na.exclude")
print(htest)
```

```
## 
## 	Wilcoxon signed rank test
## 
## data:  actual_dnds$dNdS and expected_dnds$Omega
## V = 749, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.4476 -0.3933
## sample estimates:
## (pseudo)median 
##        -0.4062
```


**Scatterplot actual vs expected dn ds together**


```r
fullDat <- data.frame(site = expected_dnds$Site, actual = actual_dnds$dNdS, 
    expected = expected_dnds$Omega)
ggplot(fullDat, aes(x = actual, y = expected)) + geom_smooth(method = lm)
```

```
## Warning: Removed 2446 rows containing missing values (stat_smooth).
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


**Scatterplot the dn-ds across the genome**


```r
fullDatBySource <- reshape2:::melt.data.frame(data = fullDat, na.rm = TRUE, 
    id.vars = "site", variable.name = "source", value.name = "dnds")
head(fullDatBySource)
```

```
##   site source     dnds
## 1    1 actual 4.213074
## 2    2 actual 0.008495
## 3    3 actual 0.013434
## 4    4 actual 0.113133
## 5    5 actual 0.206868
## 6    6 actual 0.128783
```

```r
tail(fullDatBySource)
```

```
##      site   source dnds
## 5995 2995 expected 0.55
## 5996 2996 expected 0.45
## 5997 2997 expected 0.05
## 5998 2998 expected 0.45
## 5999 2999 expected 0.15
## 6000 3000 expected 0.15
```

```r
str(fullDatBySource)
```

```
## 'data.frame':	3554 obs. of  3 variables:
##  $ site  : int  1 2 3 4 5 6 7 8 9 10 ...
##  $ source: Factor w/ 2 levels "actual","expected": 1 1 1 1 1 1 1 1 1 1 ...
##  $ dnds  : num  4.21307 0.00849 0.01343 0.11313 0.20687 ...
```

```r
summary(fullDatBySource)
```

```
##       site           source          dnds       
##  Min.   :   1   actual  : 554   Min.   :-0.072  
##  1st Qu.: 446   expected:3000   1st Qu.: 0.150  
##  Median :1224                   Median : 0.350  
##  Mean   :1313                   Mean   : 0.420  
##  3rd Qu.:2112                   3rd Qu.: 0.550  
##  Max.   :3000                   Max.   : 4.213
```

```r
ggplot(fullDatBySource[fullDatBySource$site < 500, ], aes(x = site, y = dnds, 
    color = source)) + geom_point() + xlab("Codon Site Along Genome") + ylab("Normalized dN-dS") + 
    ggtitle("dn-ds by site")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```r
dnds_cor <- cor(actual_dnds$dNdS, expected_dnds$Omega, method = "spearman", 
    use = "pairwise.complete.obs")
print(dnds_cor)
```

```
## [1] 0.02211
```


