# Generate the phase detail for submissions

The phase details is as simple as creating a lower case version of the
phase. However, when several phases are combined in `metadata_phase.txt`
(pipe separated), a new version has to be generated:

- Animal: `PASS1A-06|PASS1C-06` generates either `pass1ac-06` or
  `pass1ac-18`

- Human: the tranche information is concatenated in canonical order
  (`HUMAN-PRECOVID` first, tranches in ascending order), so that the
  same set of phases always generates the same phase details. For
  example:

  - `HUMAN-MAIN-TR01|HUMAN-MAIN-TR02|HUMAN-MAIN-TR03`:
    `human-main-tr01-tr02-tr03`

  - `HUMAN-PRECOVID|HUMAN-MAIN-TR01|HUMAN-MAIN-TR02`:
    `human-precovid-main-tr01-tr02`

Human and animal phases cannot be combined: it throws an error.

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
