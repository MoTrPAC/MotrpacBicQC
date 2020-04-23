# MotrpacBicQC

*An R package for the MoTrPAC community*

---

[![CircleCI](https://circleci.com/gh/MoTrPAC/MotrpacBicQC.svg?style=svg)](https://circleci.com/gh/MoTrPAC/MotrpacBicQC)

  
## Overview

This package provides a set of functions for Primary and Secondary 
QC/QA analysis of datasets generated and submitted by the 
Chemical Analysis Sites (CAS).

Currently includes Metabolomics datasets only.

## Installation

### R and RStudio

Download and install R and RStudio (if you don't have it already)

### This package

Open RStudio and

- Install `devtools` library

```
install.packages("devtools")
```

- Install this package from github

```
library(devtools)
devtools::install_github("MoTrPAC/MotrpacBicQC")
```

## Usage

First, load the library

```
library(MotrpacBicQC)
```

And run any of the following tests to check that the package 
is correctly installed it works. For example:

```
# Just copy and paste in the RStudio terminal

check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named")
check_metadata_samples(df = metadata_sample_named, cas = "umichigan")
check_results(r_m = results_named, m_s = metadata_sample_named, m_m = metadata_metabolites_named)
```

### How to test your datasets

Two approaches available:

#### Check full `PROCESSED_YYYYMMDD` folder

Run test on the full submission. For that, run the following command:

```
validate_metabolomics(input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD", 
                      cas = "your_site_code")
```

**cas** is one of the followings:

- "broad_met" = Broad Metabolomics
- "emory"     = Emory
- "mayo"      = Mayo Clinic
- "umichigan" = Umichigan
- "gtech"     = Georgia Tech
- "duke"      = Duke

#### Check individual files

- Check metadata metabolites:

```
# Open the metadata_metabolites file(s)

metadata_metabolites_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
metadata_metabolites_unnamed <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)

check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named")
check_metadata_metabolites(df = metadata_metabolites_unnamed, name_id = "unnamed")

```

- Check metadata samples:

```
# Open your files
metadata_sample_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
metadata_sample_unnamed <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)

check_metadata_samples(df = metadata_sample_named, cas = "your_side_id")
check_metadata_samples(df = metadata_sample_unnamed, cas = "your_side_id")
```

- Check results, which needs both both metadata metabolites and samples

```
# Open your files
metadata_metabolites_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
metadata_sample_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
results_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)

check_results(r_m = results_named, 
              m_s = metadata_sample_named, 
              m_m = metadata_metabolites_named)
```

## Help

Need help? Please, [submit an issue here](https://github.com/MoTrPAC/MotrpacBicQC/issues) 
providing as many details as possible.

## Credit

[MoTrPAC Bioinformatics Center](https://motrpac-data.org/)



