# Validate 'lc_column_id' column

This function checks the 'lc_column_id' column of a provided data frame
to ensure that it exists, contains no NA values, and contains no spaces
in its entries. It also reports the number of unique values in the
column.

## Usage

``` r
validate_lc_column_id(df, column_name, verbose = TRUE)
```

## Arguments

- df:

  A data frame that should contain the 'lc_column_id' column.

- column_name:

  The name of the column in `df` which contains the date-time values.

- verbose:

  A logical indicating whether to print informative messages. Default is
  TRUE.

## Value

An invisible NULL. The function is used mainly for its side effects
(i.e., printing validation results).

## Examples

``` r
df <- data.frame(lc_column_id = c("id1", "id2", "id3", "id1", "id 2", NA))
validate_lc_column_id(df, column_name = "lc_column_id")
#>    - (-) NA values detected in column ` lc_column_id `: FAIL
#>    - (-) Spaces detected in column `lc_column_id`: FAIL
#>    - ( ) Number of unique values in column ` lc_column_id `:  5
#>    - (!) Warning: the number of LC columns might be too high. Please, revise 
#> [1] 2
```
