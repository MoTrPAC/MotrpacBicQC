# Load metabolomics batch

Open, check, and return all metabolomics files

## Usage

``` r
load_metabolomics_batch(input_results_folder, cas, verbose = TRUE)
```

## Arguments

- input_results_folder:

  (char) Path to the PROCESSED_YYYYMMDD folder

- cas:

  (char) Chemical Analytical Site code (e.g "umichigan")

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
