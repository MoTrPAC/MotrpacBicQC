# Load Proteomics batch

Load a proteomics batch

## Usage

``` r
load_proteomics(input_results_folder, isPTM, verbose = TRUE)
```

## Arguments

- input_results_folder:

  (char) path to the PROCESSED folder to check

- isPTM:

  (logical) `TRUE` if it is Post-Translational-Modification proteomics
  assay

- verbose:

  (logical) `TRUE` (default) prints QC details.

## Value

(List) of data frames: rii,
