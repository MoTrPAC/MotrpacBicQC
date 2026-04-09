# QC Plot of proteomics reporter ion intensity data

check whether results file is following guidelines

## Usage

``` r
proteomics_plots_rii(
  all_vial_labels,
  all_samples,
  peprii,
  isPTM,
  v_m,
  out_qc_folder = NULL,
  output_prefix,
  printPDF,
  verbose
)
```

## Arguments

- all_vial_labels:

  (vector) only vial_labels

- all_samples:

  (vector) all sample ids

- peprii:

  (char) Reporter Ion Intensity df

- isPTM:

  (logical) Is a PTM assay

- v_m:

  (df) sample metadata

- out_qc_folder:

  (char) if `f_proof is TRUE`, a folder path can be provided (otherwise
  print in current working directory)

- output_prefix:

  (char) provide a prefix for the output name

- printPDF:

  (logical) if `TRUE` (default print plots to pdf)

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
