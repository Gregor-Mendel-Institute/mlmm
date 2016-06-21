# MLMM

## Introduction

This directory contains the `mlmm` package for the R programming language. It implements an efficient multi-locus mixed-model approach for genome-wide association studies in structured populations.

## Authors and license

The main authors are Vincent Segura and Bjarni J. Vilhjalmsson. The code is available under the GNU Public License (version 3 and later). See the COPYING file for usage permissions.

## Development

The content of this directory is versioned using git, the central repository being hosted on [GitHub](https://github.com/Gregor-Mendel-Institute/mlmm). Please report issues directly [online](https://github.com/Gregor-Mendel-Institute/mlmm/issues).

## Installation

For users, the easiest is to directly install the package from GitHub:
```
R> library(devtools); install_github("Gregor-Mendel-Institute/mlmm")
```

Note that this package depends on the `emma` package (not the one on CRAN, but the one from UCLA available [here](http://mouse.cs.ucla.edu/emma/)).

For developpers, when editing the content of this repo, increment the version of the package in `DESCRIPTION` and execute the following commands:
```
$ Rscript -e 'library(devtools); devtools::document()'
$ R CMD build mlmm
$ R CMD check mlmm_<version>.tar.gz
$ sudo R CMD INSTALL mlmm_<version>.tar.gz
```

More information is available in Hadley Wickham's [book](http://r-pkgs.had.co.nz/).

## Usage

Two main functions can be used to carry out GWAS with MLMM and plot the results from the analysis:

* `mlmm`, the original MLMM as described in [Segura, Vilhjálmsson et al. (Nat Gen 2012)](http://www.nature.com/ng/journal/v44/n7/full/ng.2314.html).

* `mlmm_cof`, a modified version of MLMM that allows including a fixed covariate in the association model. This could be for example a matrix of principal components scores (MLMM version of the "PK" model) or any feature that would make sense to regress out (e.g. sex).

In their current versions, the MLMM functions do not allow for missing values in the genotype matrix. Whenever possible we would suggest imputing the genotypic data prior to the analysis.

Once the package is installed, browse the vignettes:
```
R> library(mlmm)
R> browseVignettes("mlmm")
```

When used for a scientific article, don't forget to cite it:
```
R> citation("mlmm")
```

See also `citation()` for citing R itself.
