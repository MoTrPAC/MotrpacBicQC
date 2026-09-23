# Plot the percentage of NA values per column

Bar plot with the percentage of `NA` values in every column of a data
frame. The bars follow the order of the columns in the data frame. The
percentages are calculated with
[`naniar::miss_var_summary()`](https://naniar.njtierney.com/reference/miss_var_summary.html).

This function replaces the
`inspectdf::inspect_na() %>% inspectdf::show_plot()` chart previously
used in the QC plots: `inspectdf` was archived from CRAN (2026-04-10),
which made the package impossible to install from a clean library.

## Usage

``` r
plot_na_percentage(df, text_labels = TRUE)
```

## Arguments

- df:

  (data.frame) data frame to inspect

- text_labels:

  (logical) `TRUE` (default) prints the percentage of `NA` values on
  every bar

## Value

(ggplot) bar plot with the percentage of `NA` values per column

## Examples

``` r
plot_na_percentage(results_named)
```
