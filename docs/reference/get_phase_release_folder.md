# Get the phase folder of the release

The release folder structure is `PHASE/OMICS/TCODE_NAME/ASSAY/`, where
the `PHASE` folder is the phase details with two exceptions:

- `pass1c-06` releases are written to the `pass1a-06` folder

- Every `HUMAN-MAIN` release is written to the `human-main` folder, no
  matter how many tranches are combined, and whether `HUMAN-PRECOVID` is
  combined with them. For example, `human-main-tr01`,
  `human-main-tr01-tr02-tr03`, and `human-precovid-main-tr01-tr02` are
  all written to `human-main`

The full phase details (including every tranche) is still used for the
file names, so the content of the release folder identifies the combined
phases.

## Usage

``` r
get_phase_release_folder(phase_details)
```

## Arguments

- phase_details:

  (char) expected output of `generate_phase_details`

## Value

(char) name of the `PHASE` folder of the release
