# check metadata metabolites

check whether metadata_metabolites is following guidelines

## Usage

``` r
check_metadata_metabolites(
  df,
  name_id,
  return_n_issues = FALSE,
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
#>   + (+) `refmet_name` unique values: OK
#>   + Validating `refmet_name` (it might take some time)
#>       (-) `refmet_name` [`Leucine/Isoleucine`] must be modified to the RefMet Standarized name: "Leucine" (Error RN2)
#>       (-) `refmet_name` [`Oxoglutaric acid`] must be modified to the RefMet Standarized name: "2-Oxoglutaric acid" (Error RN2)
#>       (-) `refmet_name` [`Citric acid/Isocitric acid`] must be modified to the RefMet Standarized name: "Citric acid" (Error RN2)
#>       (-) Total number of missed ids on MW: 3
#>    - (-) SUMMARY: 3 `refmet_name` not found in RefMet Metabolomics Data Dictionary: FAIL
#>   + (+) {rt} all numeric: OK
#>   + (+) {mz} all numeric: OK
#>   + (+) {`neutral_mass`} all numeric values OK
#>   + (+) {formula} available: OK
```
