# MotrpacBicQC

*An R package for the MoTrPAC community*

---

<!-- badges: start -->
[![CircleCI](https://circleci.com/gh/MoTrPAC/MotrpacBicQC.svg?style=svg)](https://circleci.com/gh/MoTrPAC/MotrpacBicQC)
![R package version](https://img.shields.io/github/r-package/v/MoTrPAC/MotrpacBicQC?label=R%20package)
![Last commit](https://img.shields.io/github/last-commit/MoTrPAC/MotrpacBicQC/develop)
[![DOI](https://zenodo.org/badge/256275809.svg)](https://zenodo.org/badge/latestdoi/256275809)
<!-- badges: end -->

  
## Overview

This package provides a set of functions for *Primary* and *Secondary* 
QC/QA analysis of datasets generated by the Metabolomics, Targeted, and Untargeted Proteomics MoTrPAC
Chemical Analysis Sites (CAS). It also provides additional functions and data objects for analysis.

## Package Updates

This package is frequently updated. Please, check the file [NEWS.md](https://motrpac.github.io/MotrpacBicQC/news/index.html) to find out more about the changes affecting every version

## Installation

First, download and install R and RStudio:

- [R](https://mirror.las.iastate.edu/CRAN/) 
- [RStudio](https://rstudio.com/products/rstudio/download/) (free version)

Then, open RStudio and install the `devtools` package

```
install.packages("devtools")
```

Finally, install the `MotrpacBicQC` package

```
library(devtools)
devtools::install_github("MoTrPAC/MotrpacBicQC")
```

## Usage

The following vignettes are available:

- [Metabolomics](https://motrpac.github.io/MotrpacBicQC/articles/qc_metabolomics.html)
- [Proteomics (untargeted)](https://motrpac.github.io/MotrpacBicQC/articles/qc_proteomics.html)
- [Proteomics (olink)](https://motrpac.github.io/MotrpacBicQC/articles/qc_olink.html)
- [Clinical Chemistry](https://motrpac.github.io/MotrpacBicQC/articles/qc_lab.html)
- [Other functions](https://motrpac.github.io/MotrpacBicQC/articles/other_functions.html)

Alternatively, once the package is installed, run the following command to 
access the same documentation:

```
browseVignettes("MotrpacBicQC")
```

[Follow this link](https://motrpac.github.io/MotrpacBicQC/reference/index.html) 
for details of all the available functions


## Help

Need help? Please, [submit an issue here](https://github.com/MoTrPAC/MotrpacBicQC/issues) 
providing as many details as possible.


## Credit

[MoTrPAC Bioinformatics Center](https://motrpac-data.org/)



