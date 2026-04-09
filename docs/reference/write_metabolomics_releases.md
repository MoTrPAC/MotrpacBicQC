# Write metabolomics data release

Write out metabolomics data releases. Doesn't check whether data has
been submited according to guidelines

## Usage

``` r
write_metabolomics_releases(
  input_results_folder,
  cas,
  folder_name = "motrpac_release",
  folder_root = NULL,
  version_file = "v1.0",
  verbose = TRUE
)
```

## Arguments

- input_results_folder:

  (char) Path to the PROCESSED_YYYYMMDD folder

- cas:

  (char) Chemical Analytical Site code (e.g "umichigan")

- folder_name:

  (char) output files name. Must have a `.yaml` extension.

- folder_root:

  (char) absolute path to write the output files. Default: current
  directory

- version_file:

  (char) file version number (`v#.#`)

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

bic release folder/file structure `PHASE/OMICS/TCODE_NAME/ASSAY/` and
file names, including:

- `motrpac_YYYYMMDD_phasecode_tissuecode_omics_assay_file-details.txt`
  where files-details can be:

- `named-experimentalDetails.txt`

- `named-metadata-metabolites.txt`

- `metadata-samples.txt`

- `named-results.txt`

## Examples

``` r
if (FALSE) { # \dontrun{
write_metabolomics_releases(input_results_folder = "/path/to/PROCESSED_YYYYMMDD/", 
cas = "umichigan")
} # }
```
