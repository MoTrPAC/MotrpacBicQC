# Metabolomics data dictionary

Metabolomics data dictionary

## Usage

``` r
metabolomics_data_dictionary
```

## Format

A data frame with \> 2000 rows and 15 columns (variables):

- CURRENT_REFMET_NAME:

  Updated version of RefMet (after December 2020)

- refmet_name:

  Old version of refmet_name

- metabolite_name:

  Metabolite name provided by each laboratory

- is_standard:

  is a reference standard? 1: yes, 0: no

- super_class:

  Metabolite super class

- main_class:

  Metabolite main class

- sub_class:

  Metabolite sub class

- formula:

  Metabolite formula

- exactmass:

  Metaboliteexact mass

- pubchem_cid:

  Metabolite pubchem id

- kegg_id:

  Metabolite kegg id

- inchi_key:

  Metabolite inchi key

- lm_id:

  Metabolite lm id

- hmdb_id:

  Metabolite hmdb id

- chebi_id:

  Metabolite chebi id

## Examples

``` r
if (FALSE) { # \dontrun{
 metabolomics_data_dictionary
} # }
```
