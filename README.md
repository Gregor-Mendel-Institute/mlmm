# MLMM - An efficient multi-locus mixed-model approach
--
An efficient multi-locus mixed-model approach for genome-wide association studies in structured populations

License: GPL-3
Copyright: Vincent Segura, Bjarni J. Vilhjalmsson

## Introduction

This repository contains R functions to carry out GWAS with MLMM and plot the results from the analysis.

Two versions are currently available:

* [mlmm.r](mlmm.r), the original MLMM as described in [Segura, Vilhjálmsson et al. Nat Gen 2012](http://www.nature.com/ng/journal/v44/n7/full/ng.2314.html).
* [mlmm_cof.r](mlmm_cof.r) a modified version of MLMM that allow including a fixed covariate in the association model. This could be for example a matrix of principal components scores (MLMM version of the "PK" model) or any feature that would make sense to regress out (e.g. sex).

In their current versions, the MLMM functions do not allow for missing values in the genotype matrix. Whenever possible we would suggest imputing the genotypic data prior to the analysis.

## How to use

### MLMM.r
R script to perform mlmm

```R
#set your working directory
setwd('my_directory')
 
#load the tutorial data for carrying out mlmm analysis
load('data/example_data.Rdata')
str(example_data)
 
#load the mlmm function as well as the emma package (if it does not install with your current R version, just download and source it, as recommended on the emma website).
source('mlmm.r')
source('emma.r')
 
#perform mlmm (10 steps), it can take few minutes...
mygwas<-mlmm(Y=example_data$Y,X=example_data$X,K=example_data$K,nbchunks=2,maxsteps=10)
 
#display and plot the results
#mlmm stepwise table
mygwas$step_table
#EBIC plot
plot_step_table(mygwas,'extBIC')
#mbonf criterion plot
plot_step_table(mygwas,'maxpval')
#% variance plot
plot_step_RSS(mygwas)
#1st mlmm step plot
plot_fwd_GWAS(mygwas,1,example_data$snp_info,0.1)
#2nd mlmm step plot
plot_fwd_GWAS(mygwas,2,example_data$snp_info,0.1)
#3rd mlmm step plot
plot_fwd_GWAS(mygwas,3,example_data$snp_info,0.1)
#optimal step according to ebic plot
plot_opt_GWAS(mygwas,'extBIC',example_data$snp_info,0.1)
#optimal step according to mbonf plot
plot_opt_GWAS(mygwas,'mbonf',example_data$snp_info,0.1)
```

### MLMM_COF.r
R script to perform mlmm with covariates

```R
#set your working directory
setwd('my_directory')

#load the tutorial data for carrying out mlmm analysis
load('data/example_data_bis.Rdata')
#this object is similar to the file "example_dara.Rdata", except that it also contains PCs scores of individuals to be used as covariate in mlmm ("PK" model)
str(example_data)

#load the mlmm function as well as the emma package (if it does not install with your current R version, just download and source it, as recommended on the emma website).
source('mlmm_cof.r')
source('emma.r')

#perform mlmm (5 steps), it can take few minutes...
mygwas<-mlmm_cof(example_data$Y,example_data$X,example_data$PC[,1:10],example_data$K,10,5)

#display and plot the results
#mlmm stepwise table

mygwas$step_table

#EBIC plot
plot_step_table(mygwas,'extBIC')
#mbonf criterion plot
plot_step_table(mygwas,'maxpval')
#% variance plot
plot_step_RSS(mygwas)
#1st mlmm step plot
plot_fwd_GWAS(mygwas,1,example_data$snp_info,0.1)
#2nd mlmm step plot
plot_fwd_GWAS(mygwas,2,example_data$snp_info,0.1)
#3rd mlmm step plot
plot_fwd_GWAS(mygwas,3,example_data$snp_info,0.1)
#optimal step according to ebic plot
plot_opt_GWAS(mygwas,'extBIC',example_data$snp_info,0.1)

#optimal step according to mbonf plot
plot_opt_GWAS(mygwas,'mbonf',example_data$snp_info,0.1)

#qqplots for 5 steps and for optimal models
qqplot_fwd_GWAS(mygwas,5)
qqplot_opt_GWAS(mygwas,'extBIC')
qqplot_opt_GWAS(mygwas,'mbonf')
```
