# check metadata metabolites

check whether metadata_metabolites is following guidelines

## Usage

``` r
check_metadata_metabolites(
  df,
  name_id,
  return_n_issues = FALSE,
  refmet_validation = TRUE,
  verbose = TRUE
)
```

## Arguments

- df:

  (data.frame) metadata_metabolites

- name_id:

  (char) specify whether `named` or `unnamed` files

- return_n_issues:

  (logical) if `TRUE` returns the number of issues.

- refmet_validation:

  (logical) `TRUE` (default) validates every `refmet_name` against the
  Metabolomics Workbench API (one request per metabolite, slow). `FALSE`
  skips only the API calls: the `refmet_name` column-presence and
  uniqueness checks always run.

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
{
check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named")
}
#>   + (+) All required columns present
#>   + (+) `metabolite_name` OK
#>   + (+) `refmet_name validation`: unique values: OK
#>   + (+) `refmet_name validation`: connecting to Metabolomics Workbench for validation (slow)
#>    - (-) `refmet_name validation`: [`Oxoglutaric acid`] must be modified to the RefMet Standardized name: "2-Oxoglutaric acid" (Error RN2)
#>    - (-) `refmet_name validation`: Total number of missed ids on MW: 1
#>    - (-) `refmet_name validation`: 1 `refmet_name` not found in RefMet: FAIL
#>   + (+) {rt} all numeric: OK
#>   + (+) {mz} all numeric: OK
#>   + (+) {`neutral_mass`} all numeric values OK
#>   + (+) {formula} available: OK
```
