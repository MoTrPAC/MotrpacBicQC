---
title: "MotrpacBicQC: CHEMICAL CHEMISTRY LAB QC"
date: "2025-03-15"
output:
  rmdformats::downcute:
    code_folding: show
    self_contained: true
    thumbnails: false
    lightbox: true
pkgdown:
  as_is: true
  
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{MotrpacBicQC: OLINK QC}
  %\VignetteEncoding{UTF-8}
---


## Installation

First, download and install R and RStudio:

- [R](https://mirror.las.iastate.edu/CRAN/) 
- [RStudio](https://rstudio.com/products/rstudio/download/) (free version)

Then, open RStudio and install the `devtools` package


``` r
install.packages("devtools")
```

Finally, install the `MotrpacBicQC` package


``` r
library(devtools)
devtools::install_github("MoTrPAC/MotrpacBicQC", build_vignettes = TRUE)
```


## Usage

Load the library


``` r
library(MotrpacBicQC)
```

And run any of the following tests to check that the package 
is correctly installed and it works. For example:


``` r
# Just copy and paste in the RStudio terminal. 
check_metadata_samples_lab(df = metadata_metabolites_named)
check_metadata_analyte(df = metadata_metabolites_named)
check_results_assays(df = results_named, assay_type = "lab")
```

which should generate the following outputs:

```
   - (-) `metadata_samples`: Expected COLUMN NAMES are missed: FAIL
	 The following required columns are not present: `sample_id, sample_type, sample_order, raw_file, extraction_date, acquisition_date`
   - (-) `sample_id` column missing: FAIL
   - (-) `sample_type` column missing: FAIL
   - (-) `sample_order` column missing: FAIL
   - (-) `raw_file` column missing: FAIL
   - (-) `extraction_date` column missed: FAIL
   - (-) `acquisition_date` column missed: FAIL
   - (-) `metadata_analytes`: Expected COLUMN NAMES are missed: FAIL
	 The following required columns are not present: `analyte_name, uniprot_entry, assay_name`
   - (-) `analyte_name` column missing: FAIL
   - (-) `uniprot_entry` column missing: FAIL
   - (-) `assay_name` column missing: FAIL
   - (-) `analyte_name` column missing: FAIL
   - (-) `results` contains non-numeric columns: FAIL
		 - metabolite_name
  + ( ) Number of zeros in dataset: 14 (out of 5099 values)
  + ( ) Number of NAs in dataset: 95 (out of 5194 values)
```

## How to test your datasets

Check full `PROCESSED_YYYYMMDD` folder (recommended). The typical folder and
file structure should look like this:

```
└── HUMAN
    └── T02
        ├── LAB_CK
        │   ├── BATCH1_20221102
        │   │   ├── PROCESSED_20221102
        │   │   │   ├── metadata_analyte_named_CK_plasma.txt
        │   │   │   ├── metadata_experimentalDetails_named_duke_ClinChem.txt
        │   │   │   ├── metadata_sample_named_CK_plasma.txt
        │   │   │   └── results_CK_plasma.txt
        │   │   ├── metadata_failedsamples_20221102.txt
        │   │   └── metadata_phase.txt
        │   │   └── file_manifest_20240103.csv
```


Run test on the full submission. For that, run the following command:


``` r
n_issues <- validate_lab(input_results_folder = "/full/path/to/HUMAN/T02/LAB_CK/BATCH1_20221102/PROCESSED_20221102/", 
                         cas = "duke",
                         return_n_issues = TRUE,
                         verbose = TRUE)
```

A typical output would look like this:

```
# LAB Assay QC Report


+ Site: duke  
+ Folder: `HUMAN/T02/LAB_CK/BATCH1_20221102/PROCESSED_20221102`
+ Motrpac phase reported: HUMAN-PRECOVID (info from metadata_phase.txt available): OK                

## QC `metadata_analyte` file

  + (+) File successfully opened
  + (+) All required columns present
  + (+) `analyte_name` unique values: OK
  + (+) `uniprot_entry` unique values: OK
  + Validating `uniprot_entry` IDs with the Uniprot database. Please wait...
  + (+) All `uniprot_entry` IDs are valid: OK
  + (+) `assay_name` unique values: OK

## QC `metadata_sample` file

  + (+) File successfully opened
  + (+) All required columns present
  + (+) `sample_id` unique values: OK
  + (+) `sample_type` values are valid: OK
  + (+) `sample_order` is numeric: OK
  + (+) `raw_file` values are valid: OK
  + (+) `extraction_date`: All dates are valid.
  + (+) `acquisition_date`: All dates are valid.

## QC `results` file

  + (+) File successfully opened
  + (+) `analyte_name` unique values: OK
  + (+) All measurement columns are numeric: OK
  + ( ) Number of zeros in dataset: 0 (out of 1438 values)
  + ( ) Number of NAs in dataset: 0 (out of 1438 values)

## Cross-File Validation

  + (+) All sample IDs match between results and metadata samples: OK
  + (+) All analyte IDs match between results and metadata analytes: OK

## QC Plots

  + (p) Plot QC plots: OK

## QC `file_manifest_YYYYMMDD.csv` (required)

  + (+) `file_name, md5` columns available in manifest file
  + (+) `metadata-proteins` file included in manifest: OK
  + (+) `metadata-samples` file included in manifest: OK
  + (+) `results` file included in manifest: OK


## DMAQC validation

   + ( ) File [`metadata_failedsamples.*.txt`] not found
   + ( ) NO FAILED SAMPLES reported

TOTAL NUMBER OF ISSUES: 0
```

## Help

Additional details for each function can be found by typing, for example:


``` r
?validate_lab
```

Need extra help? Please, [submit an issue here](https://github.com/MoTrPAC/MotrpacBicQC/issues) 
providing as many details as possible.

