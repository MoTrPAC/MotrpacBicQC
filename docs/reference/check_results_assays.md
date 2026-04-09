# Check results file for assays

Checks whether the results file follows guidelines for Olink and IMM
assays.

## Usage

``` r
check_results_assays(
  df,
  assay_type = c("olink", "lab"),
  return_n_issues = FALSE,
  verbose = TRUE
)
```

## Arguments

- df:

  (data.frame) The results data frame to check.

- assay_type:

  (character) The type of assay, either `"olink"` or `"imm"`.

- return_n_issues:

  (logical) If `TRUE`, returns the number of issues identified.

- verbose:

  (logical) If `TRUE` (default), displays messages during the checking
  process.

## Value

(int) Number of issues identified if `return_n_issues` is `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
check_results_assays(df = results_df, assay_type = "olink")
check_results_assays(df = results_df, assay_type = "imm")
} # }
```
