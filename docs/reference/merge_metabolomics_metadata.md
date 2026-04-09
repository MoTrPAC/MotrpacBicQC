# Merge metabolomics metadata named and unnamed

Merge metabolomics metadata

## Usage

``` r
merge_metabolomics_metadata(m_m_n, m_m_u)
```

## Arguments

- m_m_n:

  (df) metabolomics metadata named

- m_m_u:

  (char) metabolomics metadata unnamed

## Value

(data.frame) merged metadata metabolites

## Examples

``` r
{
m_m <- merge_metabolomics_metadata(m_m_n = metadata_metabolites_named,
                                   m_m_u = metadata_metabolites_unnamed)
}
```
