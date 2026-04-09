# check results

check whether results file is following guidelines

## Usage

``` r
check_results(r_m, m_s, m_m, return_n_issues = FALSE, verbose = TRUE)
```

## Arguments

- r_m:

  (data.frame) results df

- m_s:

  (char) metadata_sample df

- m_m:

  (char) metadata_metabolites df

- return_n_issues:

  (logical) if `TRUE` returns the number of issues

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified

## Examples

``` r
{
check_results(r_m = results_named, m_s = metadata_sample_named, m_m = metadata_metabolites_named)
}
#>   + (+) All samples from `results_metabolite` are available in `metadata_sample`
#>   + (+) `metabolite_name` is identical in both [results] and `metadata_metabolites` files: OK
#>   + (+) `sample_id` columns are numeric: OK
```
