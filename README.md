---
output: 
  html_document:
    keep_md: true
---



# SV-network
Joint identification of spatially variable (SV) genes via a network-assisted Bayesian regularization approach

## Introduction
SV-network is a Bayesian regularization approach for SV gene identification with the network (dependency) structure among genes well accomodated. The confounding variations introduced by latent cellular composition within spots are also effectively corrected for higher detection accuracy.

The following R packages are required.

- ``Rcpp``
- ``RcppParallel``
- ``foreach``
- ``truncnorm``
- ``pscl``
- ``igraph``
- ``MCMCpack``
- ``Brobdingnag``

## Run SV-network on example data
An example data is utilized to introduce the usage of the work.

### Load SV-network

```r
source("R/SV-network.R")
```

### Load demo data
SV-network requires four specified inputs.

1. a raw count matrix with each row representing one spot, and each column representing one gene.

2. a matrix of spatial gene expression variation with pre-specified spatialization function, typical usage is to have the program compute the matrix based on the location matrix `loc`, specified pattern `pattern`, and quantile value `q`.

3. a cellular composition matrix with each row representing the cellular composition of one spot.

4. the input network of genes.


```r
load("data/demo_DataList.Rdata")
count <- demo_DataList$count
loc <- demo_DataList$loc
W <- demo_DataList$W
net <- demo_DataList$net
```


```r
head(count[,1:10])
```

```
##     gene1 gene2 gene3 gene4 gene5 gene6 gene7 gene8 gene9 gene10
## 1x1     0     0     0     0     3     0     0     0     0      0
## 1x2     2     5     2     0     8     0     3     0     3      1
## 1x3     0     0     1     0     0     0     0     0     0      0
## 1x4     2     0    21     0     6    26     0     3     4      8
## 1x5     0     0     4     4     6    13     0     0     0      1
## 1x6     0     0     0     4     0     1     0     0     0      0
```


```r
head(loc)
```

```
##     x y
## 1x1 1 1
## 1x2 1 2
## 1x3 1 3
## 1x4 1 4
## 1x5 1 5
## 1x6 1 6
```


```r
head(W)
```

```
##           [,1]      [,2]       [,3]        [,4]       [,5]       [,6]
## 1x1 0.43789103 0.1626681 0.03111525 0.081592530 0.03363779 0.11448268
## 1x2 0.12320358 0.0186340 0.07421685 0.163099461 0.22871728 0.05720026
## 1x3 0.02683047 0.1318463 0.04013802 0.140616508 0.37613897 0.26335371
## 1x4 0.10375804 0.2845314 0.61769139 0.008190982 0.05873937 0.16446963
## 1x5 0.20589794 0.1586651 0.31293517 0.175838933 0.17027408 0.26601696
## 1x6 0.10241893 0.1528058 0.02698858 0.430661587 0.28113009 0.21985184
```

### run the SV-network model

```r
res <- SV_network(Y=count,B=NULL,W=W,net=net,pattern="linear",q=0.5,loc=loc,a=0.05,vrbs=T)
```

```
## [1] "MCMC procedure for initialization starts!"
## [1] "10% has been done"
## [1] "20% has been done"
## [1] "30% has been done"
## [1] "40% has been done"
## [1] "50% has been done"
## [1] "60% has been done"
## [1] "70% has been done"
## [1] "80% has been done"
## [1] "90% has been done"
## [1] "100% has been done"
## [1] "MCMC procedure for initialization has been done :)"
## Time difference of 33.94594 secs
## [1] "10% has been done"
## [1] "20% has been done"
## [1] "30% has been done"
## [1] "40% has been done"
## [1] "50% has been done"
## [1] "60% has been done"
## [1] "70% has been done"
## [1] "80% has been done"
## [1] "90% has been done"
## [1] "100% has been done"
## Time difference of 1.236179 mins
```

Return the set of identified SV genes under the desired significance level `a`.

```r
res$gene_ids
```

```
##  [1]   1   2   3   4   5   6   7  10  11  12  13  14  15  16  17  18  19  20  21
## [20]  22  23  24  25  26  27  29  30  31  32  33  35  36  37  39  40  41  43  44
## [39]  45  47  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65
## [58]  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80  81  82  83  84
## [77]  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 241 431
```
