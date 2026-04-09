# Create folder

Create a directory if it doesn't exist. If no argument is provided, it
returns the current working directory

## Usage

``` r
create_folder(folder_name = NULL, verbose = FALSE)
```

## Arguments

- folder_name:

  (chr) folder name

- verbose:

  (logical) `TRUE` shows messages (default `FALSE`)

## Examples

``` r
{
create_folder(folder_name = NULL)
# Or use this one for a real folder:
# create_folder(folder_name = "testing")
}
#> [1] "/Users/davidjm/github/MoTrPAC/MotrpacBicQC/docs/reference"
```
