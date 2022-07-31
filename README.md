# mtradeR

Metagenomic TRajectory Analysis with Disease Endpoint and Risk factors.

## Introduction

We proposed a Joint model with Matching and Regularization (JMR) to detect OTUs predictive of hosts' disease status in matched sets. We also designed and implemented a simulation pipeline to generate disease outcome and temporal high-dimensional metagenomic counts in matched sets. A preprint is available on [BioRxiv.](https://www.biorxiv.org/content/10.1101/2022.04.19.488854v2)

## Installation

``` r
install.packages("devtools")
library(devtools)
install_github('qianli10000/mtradeR')
library('mtradeR')
```

## A real data example

We use the real data from a small subset of longitudinal metagenomes (n=153) in a large-cohort study of type 1 diabetes [Stewart and Ajami et al.](https://www.nature.com/articles/s41586-018-0617-x) to demonstrate the usage of JMR. We randomly select P=100 OTUs with relative abundance \>$10^{-5}$ and prevalence \>10% in this example. The OTU names are masked in this dataset.

Load real data.

``` r
data(example_data)
long<-example_data$long
logistic<-example_data$logistic
taxa_filtered<-example_data$taxa_filtered
```

Preview of real data

``` r
head(long)
  subjectID setID order sampleID sample_age outcome genotype
1    274433  1161     1       S1          4       0        1
2    274433  1161     1       S2          5       0        1
3    274433  1161     1       S3          6       0        1
4    274433  1161     1       S4          7       0        1
5    274433  1161     1       S5          8       0        1
6    893569  1161     2       S6          4       1        0

head(logistic)
  subjectID setID order outcome genotype
1    274433  1161     1       0        1
2    893569  1161     2       1        0
3    290463  1217     1       0        1
4    819136  1217     2       1        1
5    592349  1308     1       0        0
6    956178  1308     2       1        1

head(taxa_filtered,c(5,5))
          S1    S2      S3    S4      S5
OTU1 0.00000 3e-05 0.00000 0e+00 0.00001
OTU2 0.00001 1e-04 0.00034 0e+00 0.00000
OTU3 0.00015 6e-04 0.00013 2e-05 0.00018
OTU4 0.00001 6e-05 0.00000 0e+00 0.00001
OTU5 0.00001 1e-05 0.00000 1e-05 0.00001
```

Format input data.

```{r}
long_design=model.matrix(~sample_age+genotype,long)
logistic_design=model.matrix(~genotype,logistic)
long_idset=long[,1:3]
logistic_idset=logistic[,1:3]
```

An intercept test without adjusting for top-correlated taxa (JMR-NC).

``` r
JMRNC.res<- JMR(otu_tab = taxa_filtered,long_design = long_design,
               logistic_design = logistic_design,outcome = logistic$outcome, 
               long_idset = long_idset,logistic_idset = logistic_idset,
               rand.var = '(Intercept)', cov.taxa=FALSE)
```

An intercept test adjusting for top-correlated taxa.

``` r
JMR.res<- JMR(otu_tab = taxa_filtered,long_design = long_design,
             logistic_design = logistic_design,outcome = logistic$outcome, 
             long_idset = long_idset,logistic_idset = logistic_idset,
             rand.var = '(Intercept)', cov.taxa=TRUE)
```

-   **otu_tab**: a relative abundance table with rows as filtered OTUs and columns as samples.
-   **long_design**: a design matrix for longitudinal variables in the abundance-presence sub-model.
-   **logistic_design**: a design matrix for risk factors in the disease sub-model.
-   **outcome**: a vector of disease outcome per subject mapped to logistic_design.
-   **long_idset**: a dataframe of subject and set identifiers mapped to long_design, in the order of subjectID, setID, order (within-set indicator).
-   **logistic_idset**: a dataframe of subject and set identifier mapped to logistic_design, in the order of subjectID, setID, order (within-set indicator).
-   **rand.var**: the type of trajectory analysis, intercept: '(Intercept)', or slope: 'sample_age'.
-   **tune**: a scalar or vector of tuning parameter for L2 regularization. If otu_tab contains \<10 OTUs, tuning is not applicable and must be set as scalar. The default value is tune=0.15.
-   **cov.taxa**: whether to adjust for unknown dependence between taxa, default is cov.taxa=TRUE. If otu_tab contains only one OTU, cov.taxa must be FALSE.
-   **n.cores**: \# of workers registered in parallel computing. If n.cores is not specified (NULL), JMR sets n.cores=detectCores()-1.

The returned dataframe test.result shows the p-value (joint.pvalue) and BH-adjusted p-value (FDR) for taxon-specific association with disease outcome.

```{r}
tail(JMRNC.res$test.result,c(5,2))
       joint.pvalue       FDR
OTU96     0.2409968 0.3596967
OTU97     0.6072212 0.6979554
OTU98     0.6588216 0.7239798
OTU99     0.6317104 0.7120273
OTU100    0.7983006 0.8492559
```

## A simulation example

Below is an example code running the simulation pipeline and testing JMR for different types of analysis.

Load baseline parameter and generate set indicator and disease outcome, format design matrix.

```{r}
data("DM_MLE")
meta_data<-StatSim(n=100)
meta_data<-meta_data[order(meta_data$set,meta_data$id),]
subj_data<-unique(meta_data[,c('id','outcome')])
outcome<-subj_data$outcome
names(outcome)<-subj_data$id
long_design <- model.matrix(~age,meta_data)
logistic<-unique(subset(meta_data,select = -c(age,ageset.id)))
rownames(logistic)=logistic$id
logistic_design <- model.matrix(~genetic,logistic)
long_idset <- meta_data[,c('id','set','order')]
logistic_idset <- logistic[,c('id','set','order')]
```

Generate metagenomic raw counts table

```{r}
raw.counts=TaxaSim(DM_MLE,StatSim = meta_data,shift_subject = 0.9,trace =F)
rel.abun=t(t(raw.counts)/colSums(raw.counts))
```

Filter taxa by relative abundance and prevalence

```{r}
mean.rel.abun=rowMeans(rel.abun)
filter=mean.rel.abun>1e-5 & rowSums(rel.abun==0)<0.9*ncol(rel.abun)
input_tab=rel.abun[filter,]
```

Run JMR for intercept test without covariate taxa and tuning, setting shrinkage at 0.15

``` r
JMR.res=JMR(otu_tab = input_tab,long_design = long_design,
            logistic_design = logistic_design,outcome = outcome, 
            long_idset = long_idset,logistic_idset = logistic_idset,
            rand.var = '(Intercept)', tune=0.15,cov.taxa=FALSE)
```

Run JMR for slope test with covariate taxa but without tuning, setting shrinkage at 0.15

``` r
JMR.res=JMR(otu_tab = input_tab,long_design = long_design,
            logistic_design = logistic_design,outcome = outcome, 
            long_idset = long_idset,logistic_idset = logistic_idset,
            rand.var = 'age', tune=0.15,cov.taxa=TRUE)
```

Run JMR for intercept test with tuning but without selecting covariate taxa

``` r
JMR.res=JMR(otu_tab = input_tab,long_design = long_design,
            logistic_design = logistic_design,outcome = outcome, 
            long_idset = long_idset,logistic_idset = logistic_idset,
            rand.var = '(Intercept)', tune=seq(0.05,0.15,0.05),cov.taxa=FALSE)
```
