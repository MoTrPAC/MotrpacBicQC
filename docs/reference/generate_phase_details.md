# Generate the phase detail for submissions

The phase details is as simple as creating a lower case version of the
phase. However, in case of PASS1A/1C a new version has to be generated:
pass1ac-06 This function detects whether there are two phases, and if
so, generate the expected version: either pass1ac-06 or pass1ac-18

## Usage

``` r
generate_phase_details(phase_metadata, verbose = TRUE)
```

## Arguments

- phase_metadata:

  (char) expected output of `set_phase`

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(char) the expected phase_details function
