---
title: "MotrpacBicQC: OLINK QC"
date: "2024-11-17"
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
  %\usepackage[UTF-8]{inputenc}
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
check_metadata_samples_olink(df = metadata_metabolites_named)
check_metadata_proteins(df = metadata_metabolites_named)
check_results_assays(df = results_named, assay_type = "olink")
```

which should generate the following outputs:

```
   - (-) `metadata_samples`: Expected COLUMN NAMES are missed: FAIL
	 The following required columns are not present: `sample_id, sample_type, sample_order, plate_id`
   - (-) `sample_id` is missed: FAIL
   - (-) `sample_type` column missed: FAIL
   - (-) `plate_id` column missed: FAIL
   - (-) `metadata_proteins`: Expected COLUMN NAMES are missed: FAIL
	 The following required columns are not present: `olink_id, uniprot_entry, assay, missing_freq, panel_name, panel_lot_nr, normalization`
   - (-) `olink_id`: is missed: FAIL
   - (-) `uniprot_entry` column missed: FAIL
   - (-) `assay` column missed: FAIL
   - (-) {missing_freq} column missed: FAIL
   - (-) `panel_name` column missed: FAIL
   - (-) `panel_lot_nr` column missed: FAIL
   - (-) `normalization` column missed: FAIL
   - (-) `olink_id` is missed: FAIL
   - (-) `results.txt` contains non numeric columns: FAIL
		 - metabolite_name
  + ( ) Number of zeros in dataset: NA (out of 5194 values)
  + ( ) Number of NAs in dataset: 95 (out of 5194 values)
```

## How to test your datasets

Check full `RESULTS_YYYYMMDD` folder (recommended). The typical folder and
file structure should look like this:

```
HUMAN/
`-- T02
    `-- PROT_OL
        `-- BATCH1_20210825
            |-- RESULTS_20221102
            |   |-- MOTRPAC_HUMAN_T02_PROT_OL_GER_20210825_metadata-proteins.txt
            |   |-- MOTRPAC_HUMAN_T02_PROT_OL_GER_20210825_metadata-samples.txt
            |   `-- MOTRPAC_HUMAN_T02_PROT_OL_GER_20210825_results.txt
            |-- file_manifest_20240103.csv
            `-- metadata_phase.txt
```


Run test on the full submission. For that, run the following command:


``` r
n_issues <- validate_olink(input_results_folder = "/full/path/to/HUMAN/T02/PROT_OL/BATCH1_20210825/RESULTS_20221102", 
                           cas = "broad_rg",
                           return_n_issues = TRUE,
                           verbose = TRUE)
```

A typical output would look like this:

```
# OLINK QC report


+ Site: broad_rg
+ Folder: `HUMAN/T02/PROT_OL/BATCH1_20210825/RESULTS_20221102`
+ Motrpac phase reported: HUMAN-PRECOVID (info from metadata_phase.txt available): OK                                                                  

## QC `metadata_proteins`

  + (+) File successfully opened
  + (+) All required columns present
  + (+) `olink_id`: unique values: OK
   - ( ) `uniprot_entry` non-unique values detected (n duplications = 3). This is OK
		 - P01375
		 - P05231
		 - P10145
   - ( ) `assay` non-unique values detected (n duplications = 3). This is OK
	 - TNF
	 - IL6
	 - CXCL8
  + (+) {missing_freq} all numeric: OK
  + (+) {panel_name} checking available panels:
	 - Cardiometabolic
	 - Neurology
	 - Inflammation
	 - Oncology
  + (+) {panel_lot_nr} checking available panels:
	 - B04409
	 - B04410
	 - B04407
	 - B04408
  + (+) {normalization} checking available panels:
	 - Intensity

## QC `metadata-samples.txt`

  + (+) File successfully opened
  + (+) All required columns present
  + (+) `sample_id` seems OK
  + (+) `sample_type` seems OK
  + (+) `plate_id` is available: OK
  + (+) `sample_order` is numeric: OK
  + (+) All `plate_id` values have unique `sample_order` values: OK


## QC `results.txt`

  + (+) File successfully opened
  + (+) `olink_id` seems OK
  + (+) All columns (except `olink_id`) are  numeric: OK
  + ( ) Number of zeros in dataset: 55 (out of 1105472 values)
  + ( ) Number of NAs in dataset: 0 (out of 1105472 values)

## Cross File Validation

  + (+) All samples in `results.txt` are available in `metadata-samples.txt`
  + (+) All `olink_id` from `results.txt` are available in `metadata-proteins.txt`

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
?validate_olink
```

Need extra help? Please, [submit an issue here](https://github.com/MoTrPAC/MotrpacBicQC/issues) 
providing as many details as possible.

