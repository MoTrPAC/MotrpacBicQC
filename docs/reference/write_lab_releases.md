# Write LAB Data Release

Write out LAB data releases. This function doesn't check whether data
has been submitted according to guidelines; it assumes that the data has
already been validated and is ready for release.

## Usage

``` r
write_lab_releases(
  input_results_folder,
  folder_name = "motrpac_release",
  folder_root = NULL,
  version_file = "v1.0",
  verbose = TRUE
)
```

## Arguments

- input_results_folder:

  (character) Path to the `PROCESSED_YYYYMMDD` folder containing the LAB
  data.

- folder_name:

  (character) Output folder name. Default is `"motrpac_release"`.

- folder_root:

  (character) Absolute path where the output folder will be created.
  Default is the current working directory.

- version_file:

  (character) File version number (e.g., `"v1.0"`). Default is `"v1.0"`.

- verbose:

  (logical) If `TRUE` (default), displays messages during processing.

## Value

Creates the BIC release folder and file structure:
`PHASE/OMICS/TCODE_NAME/ASSAY/`, and writes out files with names:
`motrpac_phase-code_tissuecode_assay_file-details-version.txt`, where
`file-details` can be:

- `metadata-analyte`,

- `metadata-samples`,

- `results`

## Examples

``` r
if (FALSE) { # \dontrun{
write_lab_releases(
   input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD/")
} # }
```
