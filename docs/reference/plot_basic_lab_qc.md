# Generate Basic LAB QC Plots

This function creates a set of four QC plots for an LAB assay. The plots
include:

- An overall value distribution histogram with a density overlay and
  missing value percentage.

- A boxplot comparing value distributions across sample types.

- A scatter plot with a loess smooth showing the trend of values over
  injection or sample order.

- A bar plot summarizing the percentage of missing data by sample type.

## Usage

``` r
plot_basic_lab_qc(
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

  A data frame containing assay results. Expected to have one row per
  analyte, with columns such as `analyte_name` and other metadata.

- results_long:

  A long-format data frame containing sample-level data. Expected
  columns include: `sample_id`, `sample_type`, `value`, `sample_order`,
  and `sample_id_ordered`.

- out_qc_folder:

  Character. The folder where QC plot PDFs will be saved. If `NULL`, the
  function will create a default folder.

- output_prefix:

  Character. A prefix to be used in the names of output PDF files.

- printPDF:

  Logical. If `TRUE`, the function outputs the plots to PDF files.

- verbose:

  Logical. If `TRUE`, the function prints progress messages to the
  console.

## Value

Invisibly returns a
[`gridExtra::grid.arrange`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html)
object containing the arranged plots. The primary output is the saved
PDF file(s) if `printPDF = TRUE`.

## Details

The function uses `ggplot2` syntax to generate clean and informative
plots, avoiding clutter even with large numbers of samples (e.g., \>
1400 samples). It calculates an overall missing value percentage for the
measured values, and provides visual summaries of the data distribution
and missing data prevalence.
