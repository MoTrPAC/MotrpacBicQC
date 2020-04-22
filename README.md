MotrpacBicQC
===

*An R package for the MoTrPAC community*

---

[![CircleCI](https://circleci.com/gh/MoTrPAC/MotrpacBicQC.svg?style=svg&circle-token=f33574f27b84d17c137bc7630ae79fdb9e6fa301)](https://circleci.com/gh/MoTrPAC/MotrpacBicQC)

  
## Overview

This package provides a set of functions for Primary and Secondary 
QC/QA analysis of datasets generated and submitted by the 
Chemical Analysis Sites (CAS).

Currently includes Metabolomics datasets only.

## Installation

Only available from github

- Install `devtools` library

```
install.packages("devtools")
```

- Install from github

```
library(devtools)
install_github("MoTrPAC/MotrpacBicQC")
```

## Usage

Load the library

```
library(MotrpacBicQC)
```

And test some of the functions with the available test datasets. For example:

```
check_metadata_metabolites(df = metadata_metabolites_named, nameun = "named")
check_metadata_samples(df = metadata_sample_named, cas = "umichigan")
check_results(r_m = results_named, m_s = metadata_sample_named, m_m = metadata_metabolites_named)
```


## Notes for developers

Recommended steps using Rstudio (requires to set up `roxygen2` & RStudio)

- Build documentation (`command + shift + D`)
- Run the tests (`command + shift + T`)
- Check everything (`command + shift + E`), which also run tests

To load the package (`command + shift + D`)

The R CMD check report should look like this:

```
── R CMD check results ───────────────────────────────── MotrpacBicQC 0.1.0 ────
Duration: 41.6s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded
```


