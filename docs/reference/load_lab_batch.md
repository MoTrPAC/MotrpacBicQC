# Load and Process LAB Batch Data

This function loads LAB batch data from the specified input directory.
It performs quality checks on the data and loads specific files related
to LAB data, including metadata for analytes and samples, and the
results file. It also integrates validation checks and warns if there
are too many issues identified in the data.

## Usage

``` r
load_lab_batch(input_results_folder, verbose = TRUE)
```

## Arguments

- input_results_folder:

  A string representing the path to the folder containing LAB batch data
  to be loaded and processed.

- verbose:

  Logical; if `TRUE`, prints detailed messages during the loading
  process.

## Value

A list containing data frames for metadata of analytes (`m_a`), metadata
of samples (`m_s`), and LAB results (`r_o`). If certain files are not
available, the corresponding entries in the list will be `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
list_of_df <- load_lab_batch(input_results_folder = "/path/to/PROCESSED_YYYYMMDD/")
} # }
```
