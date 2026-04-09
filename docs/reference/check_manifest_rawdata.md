# check rawfiles between manifest and metabolite_sample matches

check rawfiles between manifest and metabolite_sample matches

## Usage

``` r
check_manifest_rawdata(
  input_results_folder,
  m_s_n_raw,
  return_n_issues = FALSE,
  verbose = TRUE
)
```

## Arguments

- input_results_folder:

  (char) input path folder

- m_s_n_raw:

  (list) list of raw files available in the metadata sample file

- return_n_issues:

  (logical) if `TRUE` returns the number of issues

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) number of issues identified
