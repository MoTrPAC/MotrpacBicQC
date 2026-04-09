# Validate date-time format in a data frame column

This function validates that all the values in a specified column of a
given data frame adhere to the date-time format: "MM/DD/YYYY HH:MM:SS
AM/PM". If any value does not comply with this format, it prints out
those values.

## Usage

``` r
validate_dates_times(df, column_name, verbose = TRUE)
```

## Arguments

- df:

  A data frame containing the column to be validated.

- column_name:

  The name of the column in `df` which contains the date-time values.

- verbose:

  Logical. If TRUE, messages are printed to the console.

## Value

This function returns the number of issues detected

## Examples

``` r

df <- data.frame(id = 1:6,
                 datetime = c("12/31/2023 11:59:59 PM", "01/01/2024 00:00:00 AM", 
                               "02/29/2024 12:00:00 PM", "13/01/2024 01:00:00 PM", 
                               "02/28/2024 24:00:00 PM", "02/30/2024 12:00:00 PM"))
validate_dates_times(df, "datetime", TRUE)
#>    - (-)`datetime`: Values in incorrect format: `13/01/2024 01:00:00 PM, 02/30/2024 12:00:00 PM`
#> [1] 1
```
