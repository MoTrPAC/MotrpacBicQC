# Combines all files from a Metabolomics batch

Combines all the files from an untargeted assay submission including
metadata_sample, metadata_metabolites, results for both "NAMED" and
"UNNAMED" folders

## Usage

``` r
combine_metabolomics_batch(input_results_folder, cas, verbose = TRUE)
```

## Arguments

- input_results_folder:

  (char) Path to the PROCESSED_YYYYMMDD folder

- cas:

  (char) Chemical Analytical Site code (e.g "umichigan")

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
if (FALSE) { # \dontrun{
all_datasets <- combine_metabolomics_batch(
                        input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD/",
                        cas = "umichigan")
} # }
```
