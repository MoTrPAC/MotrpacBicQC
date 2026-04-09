# Check LAB metadata samples file

Checks whether the metadata_sample file for LAB assays follows the
required guidelines.

## Usage

``` r
check_metadata_samples_lab(df, return_n_issues = FALSE, verbose = TRUE)
```

## Arguments

- df:

  (data.frame) LAB metadata samples data

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
check_metadata_samples_lab(df = metadata_sample_named_CK_plasma)
} # }
```
