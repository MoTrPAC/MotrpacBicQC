# check proteomics reported ion intensity file

check whether the proteomics rri results files is following guidelines

## Usage

``` r
check_rii_proteomics(
  df_rri,
  isPTM,
  f_proof = TRUE,
  output_prefix = "rii-file",
  out_qc_folder = NULL,
  return_n_issues = FALSE,
  printPDF = TRUE,
  verbose = TRUE
)
```

## Arguments

- df_rri:

  (data.frame) proteomics rri data frame (required)

- isPTM:

  (logical) `TRUE` if a PTM results file

- f_proof:

  (logical) `TRUE` (default) to print out distribution and NA plots

- output_prefix:

  (char) if `f_proof = TRUE`, provide a prefix for the output name

- out_qc_folder:

  (char) if `f_proof is TRUE`, a folder path can be provided (otherwise
  print in current working directory)

- return_n_issues:

  (logical) if `TRUE` returns the number of issues.

- printPDF:

  (logical) if `TRUE` (default print plots to pdf)

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
{
test <- check_rii_proteomics(df_rri = metadata_metabolites_named,
isPTM =  TRUE, return_n_issues = TRUE, verbose = FALSE)
# "test" should be NULL
}
```
