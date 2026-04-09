# Get and Validate the Entire RefMet Database from Metabolomics Workbench

This function fetches and validates the Metabolomics Data Dictionary
from the Metabolomics Workbench. It provides options to remove
duplicates.

## Usage

``` r
get_and_validate_mdd(remove_duplications = FALSE, verbose = TRUE)
```

## Arguments

- remove_duplications:

  Logical; if `TRUE`, removes duplicate entries based on the
  `refmet_name` column.

- verbose:

  Logical; if `TRUE` (default), displays progress messages and warnings
  during the function execution.

## Value

Returns a data frame with the following columns:

- `refmet_name`:

  Character; the name standarized refmet name

- `pubchem_cid`:

  Character; the PubChem compound ID.

- `lm_id`:

  Character; the LIPID MAPS ID.

- `inchi_key`:

  Character; the International Chemical Identifier Key.

- `exactmass`:

  Numeric; the exact mass of the metabolite.

- `formula`:

  Character; the chemical formula of the metabolite.

- `super_class`:

  Character; the superclass category of the metabolite.

- `main_class`:

  Character; the main class category of the metabolite.

- `sub_class`:

  Character; the subclass category of the metabolite.

- `hmdb_id`:

  Character; the Human Metabolome Database ID.

- `kegg_id`:

  Character; the Kyoto Encyclopedia of Genes and Genomes ID.

Each row of the data frame represents a unique metabolite entry from the
Metabolomics Workbench Data Dictionary.

## Details

This function downloads the entire RefMet database from the Metabolomics
Workbench using their REST API. The data is initially fetched in JSON
format and then converted to a data frame. The function checks for the
presence of a 'name' column in the data frame, renaming it to
'refmet_name' for consistency. It also provides an option to remove
duplicate entries based on the 'refmet_name' column. If duplicates are
found and `remove_duplications` is `FALSE`, the function will list the
duplicated IDs but will not remove them. This can be helpful for
reviewing the data quality and consistency.

## Examples

``` r
if (FALSE) { # \dontrun{
  refmet <- get_and_validate_mdd(remove_duplications = TRUE, verbose = TRUE)
  head(refmet)
} # }
```
