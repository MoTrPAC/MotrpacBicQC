# validate a phase with two phases (pass1a and 1c)

This function only works to validate two phases reported as for example,
'PASS1A-06\|PASS1C-06' using the separator '\|'

## Usage

``` r
validate_two_phases(phase_details, verbose = TRUE)
```

## Arguments

- phase_details:

  (char) expected output of `set_phase`

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(char) the expected phase_details function
