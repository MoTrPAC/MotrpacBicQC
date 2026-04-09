# Write olink data release

Write out olink data releases. Doesn't check whether data has been
submited according to guidelines

## Usage

``` r
write_olink_releases(
  input_results_folder,
  folder_name = "motrpac_release",
  folder_root = NULL,
  version_file = "v1.0",
  verbose = TRUE
)
```

## Arguments

- input_results_folder:

  (char) Path to the RESULTS_YYYYMMDD folder

- folder_name:

  (char) output folder name.

- folder_root:

  (char) absolute path to write the output folder. Default: current
  directory

- version_file:

  (char) file version number (v#.#)

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

bic release folder/file structure `PHASE/OMICS/TCODE_NAME/ASSAY/` and
file names, including:

- `motrpac_phase-code_tissuecode_assay_file-details-version.txt` where
  files-details can be:

- `metadata-proteins`,

- `metadata-samples`,

- `results`

## Examples

``` r
if (FALSE) { # \dontrun{
write_olink_releases(
   input_results_folder = "/full/path/to/RESULTS_YYYYMMDD/")
} # }
```
