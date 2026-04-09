# Validate Olink Data

This function validates Olink data based on a set of criteria, including
folder structure, metadata, and file contents. It supports optional
functionalities like creating a PDF report and DMAQC validation (data
only available at the BIC)

## Usage

``` r
validate_olink(
  input_results_folder,
  cas,
  return_n_issues = FALSE,
  full_report = FALSE,
  f_proof = FALSE,
  printPDF = FALSE,
  out_qc_folder = NULL,
  dmaqc_shipping_info = NULL,
  dmaqc_phase2validate = FALSE,
  validate_uniprot = FALSE,
  verbose = TRUE
)
```

## Arguments

- input_results_folder:

  A string representing the path to the folder containing Olink results
  to be validated.

- cas:

  A character string indicating the CAS number.

- return_n_issues:

  Logical; if `TRUE`, the function returns the number of detected
  issues.

- full_report:

  Logical; if `TRUE`, generates a full report of the validation process.

- f_proof:

  Logical; if `TRUE`, generates proof plots for data validation.

- printPDF:

  Logical; if `TRUE` and `f_proof` is `TRUE`, saves the plots to a PDF
  file, and in such case, then provide the desired path to output the
  PDF file in the argument `out_qc_folder`

- out_qc_folder:

  Optional; a string specifying the path to the folder where output PDF
  should be saved (only if `printPDF = TRUE`). Default: current working
  directory

- dmaqc_shipping_info:

  (char) File path to the DMAQC file. Only the BIC can use this argument

- dmaqc_phase2validate:

  (char) Provide phase to validate. This argument is not required since
  it should be extracted from the input folder or from the new required
  file `metadata_phase.txt`. Please, ignore. However, if this argument
  is provided, it will take priority and this will be the phase.

- validate_uniprot:

  Logical; if `TRUE`, validates against the UniProt database.

- verbose:

  Logical; if `TRUE`, prints detailed messages during validation.

## Value

Depending on the settings, this function may return the number of issues
found, generate reports or plots, or simply perform the validation
without returning anything.

## Examples

``` r
if (FALSE) { # \dontrun{
validate_olink("/path/to/results", cas = "broad_rg", return_n_issues = TRUE)
} # }
```
