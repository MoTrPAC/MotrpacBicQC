# filter required columns only

it returns a data frame with only the required columns for metabolomics
and proteomics

## Usage

``` r
filter_required_columns(
  df,
  type = c("m_m", "m_s", "v_m", "olproteins", "olsamples", "labanalytes", "labsamples"),
  name_id = NULL,
  verbose = TRUE
)
```

## Arguments

- df:

  (data.frame) metadata_metabolites

- type:

  (char) Type of file to filter columns:

  - `m_m`: metadata metabolites

  - `m_s`: metadata samples

  - `v_m`: proteomics vial_metadata

  - `olproteins`: olink metadata proteins

  - `olsamples`: olink metadata samples

- name_id:

  (char) specify whether `named` or `unnamed` files

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(data.frame) filtered data frame with only the required columns

## Examples

``` r
{
df_filtered <- filter_required_columns(df = metadata_metabolites_named, name_id = "named")
}
#>   + (+) All required columns present
```
