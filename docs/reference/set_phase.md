# Set the phase to be validated.

A group might choose to combine different phases, due to the
complications associated with PASS1A/1C, or to release several human
tranches together. If they choose to combine phases, the CAS must
provide a new file `metadata_phase.txt` with a single line, as for
example: `PASS1A-06|PASS1C-06` or `HUMAN-MAIN-TR01|HUMAN-MAIN-TR02`.
Human and animal phases cannot be combined. This function checks if the
file is available, and set that phase as the phases to validate. In
summary, the order of preference is:

1.  function's argument: dmaqc_phase2validate (if provided in the
    validation functions)

2.  `metadata_phase.txt` file if available in the batch folder.

3.  Phase in folder structure

## Usage

``` r
set_phase(input_results_folder, dmaqc_phase2validate, verbose = TRUE)
```

## Arguments

- input_results_folder:

  (char) path to the PROCESSED/RESULTS folder to check

- dmaqc_phase2validate:

  (data.frame) dmaqc shipping information

- verbose:

  (logical) `TRUE` (default) shows messages

## Value

(int) the phase to be validated.
