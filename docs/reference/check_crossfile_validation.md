# Cross-validate assay files (Olink and LAB)

Checks that values from the results file are available in both the
metadata analyte/protein and metadata sample files.

## Usage

``` r
check_crossfile_validation(
  r_o,
  m_s,
  m_p,
  assay_type = c("olink", "lab"),
  return_n_issues = FALSE,
  verbose = TRUE
)
```

## Arguments

- r_o:

  (data.frame) Results data frame.

- m_s:

  (data.frame) Metadata sample data frame.

- m_p:

  (data.frame) Metadata analyte/protein data frame.

- assay_type:

  (character) The type of assay, either `"olink"` or `"lab"`.

- return_n_issues:

  (logical) If `TRUE`, returns the number of issues.

- verbose:

  (logical) If `TRUE` (default), displays messages during the checking
  process.

## Value

(int) Number of issues identified if `return_n_issues` is `TRUE`.

## Examples

``` r
{
if (FALSE) { # \dontrun{
# For Olink data
check_crossfile_validation(r_o = results_olink,
                           m_s = metadata_samples_olink,
                           m_p = metadata_proteins_olink,
                           assay_type = "olink")

# For LAB data
check_crossfile_validation(r_o = results_lab,
                           m_s = metadata_samples_lab,
                           m_p = metadata_analyte_lab,
                           assay_type = "lab")
} # }
}
```
