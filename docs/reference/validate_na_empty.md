# Validate Column for NA and Empty Values

This function checks if a specified column in a data frame contains
either NA or empty values.

## Usage

``` r
validate_na_empty(df, col_name, verbose = TRUE)
```

## Arguments

- df:

  A data frame.

- col_name:

  A character string specifying the name of the column to check.

- verbose:

  A logical indicating whether to print informative messages. Default is
  TRUE.

## Value

Number of issues

## Examples

``` r
df <- data.frame(A = c("a", "", NA, "d"), B = 1:4)
validate_na_empty(df, "A")
#>    - (-) NA values detected in column ` A `: FAIL
#>    - (-) Empty values detected in column ` A `: FAIL
#> [1] 2
```
