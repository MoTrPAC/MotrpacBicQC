# Validate vial labels from DMAQC

validate vial label from DMAQC

## Usage

``` r
check_viallabel_dmaqc(
  vl_submitted,
  dmaqc_shipping_info,
  tissue_code,
  cas,
  phase,
  failed_samples,
  out_qc_folder = NULL,
  outfile_missed_viallabels,
  return_n_issues = FALSE,
  verbose = TRUE
)
```

## Arguments

- vl_submitted:

  (vector) results df

- dmaqc_shipping_info:

  (data.frame) dmaqc shipping information

- tissue_code:

  (char) tissue code

- cas:

  (char) CAS code

- phase:

  (char) phase code

- failed_samples:

  (char) metadata_metabolites df

- out_qc_folder:

  (char) output qc folder (it creates the folder if it doesn't exist)

- outfile_missed_viallabels:

  (char) file name for missed vial labels

- return_n_issues:

  (logical) if `TRUE` returns the number of issues

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified
