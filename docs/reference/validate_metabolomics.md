# Validate a Metabolomics submission

Validate a Metabolomics submission

## Usage

``` r
validate_metabolomics(
  input_results_folder,
  cas,
  return_n_issues = FALSE,
  full_report = FALSE,
  f_proof = FALSE,
  printPDF = TRUE,
  out_qc_folder = NULL,
  dmaqc_shipping_info = NULL,
  dmaqc_phase2validate = FALSE,
  refmet_validation = TRUE,
  verbose = TRUE
)
```

## Arguments

- input_results_folder:

  (char) path to the PROCESSED folder to check

- cas:

  (char) CAS code

- return_n_issues:

  (logical) if `TRUE` returns the number of issues

- full_report:

  (logical) if `FALSE` (default) it returns only the total number of
  issues. If `TRUE` returns the details of the number of issues (by
  group of files, e.g., results, metadata_metabolites, etc)

- f_proof:

  (logical) generate charts including:

  - Sample intensity distribution

  - Unique ID counts

  - NA values

- printPDF:

  (logical) if `f_proof = TRUE`, then ` printPDF = TRUE` (default)
  prints plots to pdf file and argument `out_qc_folder` should be
  provided.

- out_qc_folder:

  (char) Folder to save the pdfs (if `printPDF = TRUE`). It will create
  the folder if it doesn't exist. If this argument is not provided, the
  files will be written to the working directory

- dmaqc_shipping_info:

  (char) File path to the DMAQC file. Only the BIC can use this argument

- dmaqc_phase2validate:

  (char) Provide phase to validate. This argument is not required since
  it should be extracted from the input folder or from the new required
  file `metadata_phase.txt`. Please, ignore. However, if this argument
  is provided, it will take priority and this will be the phase.
  `metadata_phase.txt` will be ignored). Examples

  - Folder with `PASS1A-06`: type either `PASS1A-06` or leave it `NULL`

  - Both `PASS1A-06` and `PASS1C-06`: type `PASS1A-06|PASS1C-06`

  - Only `PASS1C-06`: type `PASS1C-06`

- refmet_validation:

  (logical) `TRUE` (default) validates every `refmet_name` of the named
  `metadata_metabolites` file against the Metabolomics Workbench API
  (one request per metabolite, slow). `FALSE` skips only the API calls:
  the `refmet_name` column-presence and uniqueness checks always run.
  This is the QC entry point, so validation is on by default here.

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(data.frame) Summary of issues
