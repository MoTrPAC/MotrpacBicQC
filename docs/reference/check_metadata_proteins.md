# Check olink metadata proteins file

check whether metadata-proteins.txt file is following guidelines

## Usage

``` r
check_metadata_proteins(
  df,
  return_n_issues = FALSE,
  validate_uniprot = FALSE,
  verbose = TRUE
)
```

## Arguments

- df:

  (data.frame) metadata_proteins df

- return_n_issues:

  (logical) if `TRUE` returns the number of issues.

- validate_uniprot:

  (logical) if `TRUE`, check if all uniprot ids are valid connecting to
  Uniprot. Note: depending on the number of ids, it might take several
  finish to complete the validation

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
{
if (FALSE) { # \dontrun{
check_metadata_proteins(df = metadata_proteins)
} # }
}
```
