# Plot Basic OLINK QC charts

Plot value distributions, number of unique ids, NA values

## Usage

``` r
plot_basic_olink_qc(
  results,
  results_long,
  out_qc_folder = NULL,
  output_prefix,
  printPDF = TRUE,
  verbose = TRUE
)
```

## Arguments

- results:

  (df) olink results merged with proteins metadata

- results_long:

  (df) olink results in long format and merged with sample metadata

- out_qc_folder:

  (char) output qc folder (it creates the folder if it doesn't exist)

- output_prefix:

  (char) prefix for the file name output (pdf file)

- printPDF:

  (logical) `TRUE` (default) prints pdf file

- verbose:

  (logical) `TRUE` (default) shows messages
