# Validate YYYY-MM-DD Dates in a Data Frame

Validate Dates in a Specified Column of a Data Frame

## Usage

``` r
validate_yyyymmdd_dates(df, date_column, verbose = TRUE)
```

## Arguments

- df:

  A data frame that contains the date information to be validated.

- date_column:

  A character string specifying the name of the column in `df` that
  contains the dates to be validated.

- verbose:

  A logical value indicating whether or not to print messages (default:
  `TRUE`).

## Value

number of issues found

## Details

This function checks for the validity of dates in a specified column of
a given data frame. Valid dates are in the format YYYY-MM-DD, with year
values between 2018 and 2026, month values between 1 and 12, and day
values between 1 and 31. The function prints a list of invalid dates and
a success message if all dates are valid.

## Examples

``` r
df <- data.frame(
  extraction_date = c("2022-01-31", "2023-12-01", "2025-11-30"),
  other_column = 1:3
)
ic <- validate_yyyymmdd_dates(df, "extraction_date")
#>   + (+) `extraction_date`: All dates are valid.
```
