# Validate a Proteomics submissions

Validate a Proteomics submission

## Usage

``` r
validate_proteomics(
  input_results_folder,
  isPTM,
  cas,
  dmaqc_shipping_info = NULL,
  dmaqc_phase2validate = FALSE,
  f_proof = FALSE,
  out_qc_folder = NULL,
  return_n_issues = TRUE,
  full_report = FALSE,
  printPDF = TRUE,
  verbose = TRUE,
  check_only_results = FALSE
)
```

## Arguments

- input_results_folder:

  (char) path to the PROCESSED folder to check

- isPTM:

  (logical) `TRUE` if it is Post-Translational-Modification proteomics
  assay

- cas:

  (char) CAS code

- dmaqc_shipping_info:

  (char) DMAQC file with shipping information. If not provided
  (default), then DMQAC validation is not peformed (only done at the
  BIC)

- dmaqc_phase2validate:

  (char) Provide phase to validate. This argument is not required since
  it should be extracted from the input folder or from the new required
  file `metadata_phase.txt`. However, if this argument is provided, it
  will take priority (and the phase from the input folder and the
  `metadata_phase.txt` will be ignored). Examples

  - Folder with `PASS1A-06`: type either `PASS1A-06` or leave it `NULL`

  - Both `PASS1A-06` and `PASS1C-06`: type `PASS1A-06|PASS1C-06`

  - Only `PASS1C-06`: type `PASS1C-06`

- f_proof:

  (char) print out pdf with charts including:

  - Reported Ion Intensity boxplot distribution and percentage of NA
    values per sample

  - Ratio: ratio boxplot distribution and percentage of NA values per
    samples

- out_qc_folder:

  (char) if f_proof is TRUE, then a folder must be provided

- return_n_issues:

  (logical) if `TRUE` (default) returns the number of issues

- full_report:

  (logical) if `FALSE` (default) it returns only the total number of
  issues. If `TRUE` returns the details of the number of issues (by
  group of files, e.g., results, metadata_metabolites, etc)

- printPDF:

  (logical) if `TRUE` (default print plots to pdf)

- verbose:

  (logical) `TRUE` (default) prints QC details.

- check_only_results:

  (logical) if `TRUE`, only validates results (default `FALSE`)

## Value

(data.frame) Summary of issues
