# Validate UniProt IDs

This function checks if a given vector of IDs are valid UniProt IDs. It
removes NA values, empty strings, and white spaces before validation.
Each remaining ID is then checked against the UniProt database. The
function outputs a boolean value indicating whether all IDs are valid.
It also prints messages for any ID that fails the validation.

## Usage

``` r
validate_uniprot_ids_with_uniprot(ids)
```

## Arguments

- ids:

  A character vector of potential UniProt IDs.

## Value

A boolean value. TRUE if all non-NA, non-empty, and non-whitespace IDs
in the input vector are valid UniProt IDs; FALSE otherwise.

## Note

This function requires an internet connection to access the UniProt
database. It uses the 'httr' package for HTTP requests. The function can
be slow for large datasets due to multiple web requests. Also, be aware
of potential rate limits or access restrictions on the UniProt API.

## Examples

``` r
# VALID
ids1 <- c("P12345", "Q67890", NA, "", " ")
if(validate_uniprot_ids_with_uniprot(ids1)) print("Valid UNIPROT IDs")
#> [1] "Valid UNIPROT IDs"
ids2 <- c("P12345", "Q67890", NA, "", " ", "iamwrong")
if(!validate_uniprot_ids_with_uniprot(ids1)) print("Invalid UNIPROT IDs")
```
