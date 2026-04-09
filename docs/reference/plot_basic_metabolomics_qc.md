# Plot Basic Metabolomics QC charts

Plot intensity distributions, number of unique ids, NA values and RT/MZ
density maps (only untargeted)

## Usage

``` r
plot_basic_metabolomics_qc(
  results,
  results_long,
  metametab = NULL,
  out_qc_folder = NULL,
  output_prefix,
  printPDF = TRUE,
  untargeted = TRUE,
  verbose = TRUE
)
```

## Arguments

- results:

  (df) metabolomics results (both named and unnamed already merged, if
  untargeted)

- results_long:

  (df) metabolomics results (both named and unnamed already merged, if
  untargeted, in long format)

- metametab:

  (df) metadata samples metadata (only untargeted sites)

- out_qc_folder:

  (char) output qc folder (it creates the folder if it doesn't exist)

- output_prefix:

  (char) prefix for the file name output (pdf file)

- printPDF:

  (logical) `TRUE` (default) prints pdf file

- untargeted:

  (logical) `TRUE` if the dataset is untargeted (named + unnamed
  metabolites)

- verbose:

  (logical) `TRUE` (default) shows messages
