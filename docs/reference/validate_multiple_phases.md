# validate several phases reported in `metadata_phase.txt`

Validate a pipe separated combination of phases, dispatching to the
right set of rules:

- Animal: only the `PASS1A-##|PASS1C-##` pair is supported (see
  [`validate_two_phases()`](https://motrpac.github.io/MotrpacBicQC/reference/validate_two_phases.md))

- Human: every phase must be either `HUMAN-PRECOVID` or
  `HUMAN-MAIN-TR##`, and no phase can be reported twice. Any number of
  phases is accepted.

Human and animal phases cannot be combined: it throws an error.

## Usage

``` r
validate_multiple_phases(phase_details, verbose = TRUE)
```

## Arguments

- phase_details:

  (char) expected output of `set_phase`

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(char) validation message (only if `verbose = TRUE`)
