# check proteomics vial metadata file

check whether the proteomics rri results files is following guidelines

## Usage

``` r
check_vial_metadata_proteomics(df_vm, return_n_issues = FALSE, verbose = TRUE)
```

## Arguments

- df_vm:

  (data.frame) proteomics vial_label data frame (required)

- return_n_issues:

  (logical) if `TRUE` returns the number of issues.

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
{
test <- check_vial_metadata_proteomics(df_vm = metadata_metabolites_named,
return_n_issues = TRUE, verbose = FALSE)
# "test" should be NULL
}
```
