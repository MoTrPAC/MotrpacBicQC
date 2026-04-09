# Merge all metabolomics files

Merge all metabolomics datasets, including "results" and "metadata"
files, for both targeted and untargeted datasets

## Usage

``` r
merge_all_metabolomics(m_m_n, m_m_u = NULL, m_s_n, r_m_n, r_m_u = NULL, phase)
```

## Arguments

- m_m_n:

  (metabolomics metadata named)

- m_m_u:

  (metabolomics metadata unnamed)

- m_s_n:

  (metabolomics sample named)

- r_m_n:

  (results named)

- r_m_u:

  (results unnamed)

- phase:

  (MoTrPAC Animal phase. Eg. PASS1A-06)

## Value

(data.frame) Merged data frame long format

## Examples

``` r
plasma.untargeted.merged <- merge_all_metabolomics(
       m_m_n = metadata_metabolites_named,
       m_m_u = metadata_metabolites_unnamed,
       m_s_n = metadata_sample_named,
       r_m_n = results_named,
       r_m_u = results_unnamed,
       phase = "PASS1A-06")
```
