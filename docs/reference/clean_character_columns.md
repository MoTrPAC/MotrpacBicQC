# Clean character columns

Removes leading and trailing whitespace from all character columns in a
data frame using `trimws`.

## Usage

``` r
clean_character_columns(df)
```

## Arguments

- df:

  A data frame or tibble containing character columns to be cleaned

## Value

A data frame of the same structure as the input, but with whitespace
removed from all character columns

## Details

The function uses
[`dplyr::across`](https://dplyr.tidyverse.org/reference/across.html) to
apply `trimws` to all columns that are of type character. The original
data frame structure is preserved, only the content of character columns
is modified.

## Examples

``` r
df <- data.frame(
  name = c(" John ", "Jane  ", "  Bob"),
  age = c(25, 30, 35),
  city = c("New York ", " London", " Paris ")
)
clean_df <- clean_character_columns(df)
```
