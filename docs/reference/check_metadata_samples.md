# check metadata samples

check whether metadata_sample is following guidelines

## Usage

``` r
check_metadata_samples(
  df,
  cas,
  assay = "OTHER",
  return_n_issues = FALSE,
  verbose = TRUE
)
```

## Arguments

- df:

  (data.frame) metadata_metabolites

- cas:

  (char) CAS site code

- assay:

  (char) Assay code. Options:

  - `CONV`: for Duke's conventional metabolites

  - `OTHER` (default): all the other assays.

- return_n_issues:

  (logical) if `TRUE` returns the number of issues.

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
{
check_metadata_samples(df = metadata_sample_named, cas = "umichigan")
}
#>    - (-) `metadata_sample`: recently required COLUMN NAMES are missed: Adding with NA values: FAIL
#>   + (+) All required columns present
#>   + (+) `sample_id` seems OK
#>   + (+) `sample_type` seems OK
#>   + (+) `sample_order` is numeric
#>   + (+) `sample_order` unique values OK
#>   + (+) `raw_file` unique values: OK
#>    - (-) `extraction_date` has NA values: FAIL
#>    - (-) `acquisition_date` has NA values: FAIL
#>    - (-) NA values detected in column ` lc_column_id `: FAIL
#>    - ( ) Number of unique values in column ` lc_column_id `:  1
```
