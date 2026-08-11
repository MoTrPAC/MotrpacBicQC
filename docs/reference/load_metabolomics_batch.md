# Load metabolomics batch

Open, check, and return all metabolomics files

## Usage

``` r
load_metabolomics_batch(
  input_results_folder,
  cas,
  refmet_validation = FALSE,
  verbose = TRUE
)
```

## Arguments

- input_results_folder:

  (char) Path to the PROCESSED_YYYYMMDD folder

- cas:

  (char) Chemical Analytical Site code (e.g "umichigan")

- refmet_validation:

  (logical) `FALSE` (default) skips the `refmet_name` validation against
  the Metabolomics Workbench API (one request per metabolite, slow)
  while running every other check. Set to `TRUE` to also validate the
  `refmet_name` ids. Full `refmet_name` validation belongs to
  [`validate_metabolomics()`](https://motrpac.github.io/MotrpacBicQC/reference/validate_metabolomics.md),
  which should be run before loading a batch.

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(list of data.frames) List of all the data frames

## Examples

``` r
if (FALSE) { # \dontrun{
here <- load_metabolomics_batch(input_results_folder = "/path/to/PROCESSED_YYYYMMDD/",
                                cas = "cassite")
} # }
```
