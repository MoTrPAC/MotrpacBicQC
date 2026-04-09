# extract BATCH_YYYYMMDD folder

extract BATCH_YYYYMMDD folder from input folder path. Expects the format
`BATCH#_YYYYMMDD` (batch number required). As a legacy exception,
`PROT_AC` paths also accept `BATCH_YYYYMMDD` (no batch number) to
support the historical folder
`broad/PASS1A-06/T58/PROT_AC/BATCH_20190828/`.

## Usage

``` r
validate_batch(input_results_folder)
```

## Arguments

- input_results_folder:

  (char) input_results_folder path

## Value

(vector) BATCH_YYYYMMDD folder name
