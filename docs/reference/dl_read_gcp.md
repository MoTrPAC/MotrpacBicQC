# Download and Read File from Google Cloud Storage

This function downloads a file from Google Cloud Storage (GCS) to a
local directory and reads it into R as a data frame. It uses the
`gsutil` command-line tool to handle the file download.

## Usage

``` r
dl_read_gcp(
  path,
  sep = "\t",
  header = TRUE,
  tmpdir,
  gsutil_path = "gsutil",
  check_first = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- path:

  Character. The path to the file in GCS, e.g.,
  `gs://bucket-name/file-name.csv`.

- sep:

  Character. The field separator character. Default is `\t`.

- header:

  Logical. Whether the file contains the names of the variables as its
  first line. Default is TRUE.

- tmpdir:

  Character. The local directory to which the file will be downloaded.

- gsutil_path:

  Character. The path to the `gsutil` command-line tool. Default is
  "gsutil". Now it is also supporting `gcloud` command. The full path
  should be the same as for gsutil

- check_first:

  Logical. Whether to check if the file already exists locally before
  downloading. Default is TRUE.

- verbose:

  Logical. If TRUE, prints messages about the download process. Default
  is FALSE.

- ...:

  Additional arguments passed to
  [`readr::read_delim`](https://readr.tidyverse.org/reference/read_delim.html).

## Value

A data frame containing the contents of the downloaded file.

## Details

This function first checks if the specified file exists in GCS. If the
file exists, it downloads the file to the specified local directory
(`tmpdir`). If the local directory does not exist, it will be created.
The function handles spaces in directory paths by quoting them
appropriately. If the file is successfully downloaded, it is read into R
using
[`readr::read_delim`](https://readr.tidyverse.org/reference/read_delim.html).

If the `check_first` argument is set to TRUE, the function will first
check if the file already exists locally to avoid redundant downloads.
If the file is already present locally, it will not be downloaded again.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- dl_read_gcp(
  path = "gs://bucket-name/file-name.csv",
  sep = ",",
  header = TRUE,
  tmpdir = "/local/path",
  gsutil_path = "gsutil",
  check_first = TRUE,
  verbose = TRUE
)
} # }
```
