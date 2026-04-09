# Check analyte metadata file

Checks whether the metadata_analyte file follows the required
guidelines.

## Usage

``` r
check_metadata_analyte(
  df,
  return_n_issues = FALSE,
  validate_uniprot = FALSE,
  verbose = TRUE
)
```

## Arguments

- df:

  (data.frame) The metadata_analyte data frame to check.

- return_n_issues:

  (logical) If `TRUE`, returns the number of issues found.

- validate_uniprot:

  (logical) If `TRUE`, checks if all Uniprot IDs are valid by connecting
  to the Uniprot database. Note: This may take several minutes depending
  on the number of IDs.

- verbose:

  (logical) If `TRUE` (default), displays messages during the checking
  process.

## Value

(int) The number of issues identified if `return_n_issues` is `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
check_metadata_analyte(df = metadata_analyte)
} # }
```
