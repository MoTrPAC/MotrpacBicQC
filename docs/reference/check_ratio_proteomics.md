# check proteomics ratio file

check whether the proteomics ratio results files is following guidelines

## Usage

``` r
check_ratio_proteomics(
  df_ratio,
  isPTM,
  f_proof = TRUE,
  output_prefix = "ratio-file",
  out_qc_folder = NULL,
  printPDF = TRUE,
  return_n_issues = FALSE,
  verbose = TRUE
)
```

## Arguments

- df_ratio:

  (data.frame) proteomics ratio results data frame (required)

- isPTM:

  (logical) `TRUE` if a PTM results file

- f_proof:

  (logical) `TRUE` (default) to print out distribution and NA plots

- output_prefix:

  (char) if `f_proof = TRUE`, provide a prefix for the output name

- out_qc_folder:

  (char) if `f_proof is TRUE`, a folder path can be provided (otherwise
  print in current working directory)

- printPDF:

  (logical) if `TRUE` (default print plots to pdf)

- return_n_issues:

  (logical) if `TRUE` returns the number of issues.

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
{
test <- check_ratio_proteomics(df_ratio = metadata_metabolites_named,
isPTM =  TRUE, return_n_issues = TRUE, verbose = FALSE)
# "test" should be NULL
}
```
