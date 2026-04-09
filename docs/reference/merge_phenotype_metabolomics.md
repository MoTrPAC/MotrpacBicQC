# Merge phenotypic and metabolics results

Merge phenotypic data (phenotypes_pass1a06_short) and metabolomics
merged results and metadata

## Usage

``` r
merge_phenotype_metabolomics(df_long)
```

## Arguments

- df_long:

  (data.frame) Long format of a metabolomics merged results

## Value

(data.frame) Merged file, including the following columns:

- sample_id:

  Sample Id, including vial_label and site specific QC ids

- sample_type:

  Metabolomics sample types. Check metabolomics data transfer guidelines

- sample_order:

  Order of injection on Mass Spec

- metabolite_name:

  Given name by every lab

- refmet_name:

  Map of the metabolite name to the Metabolomics RefMet database

- mz:

  mass over charge

- rt:

  retention time

- formula:

  chemical formula

- neutral_mass:

  neutral mass

- id:

  type of metabolite identification: "named", "unnamed"

- metabolite:

  Merge "refmet" for "named" metabolites and "metabolite_name" for
  "unnamed" metabolites

- quantification:

  Untargeted: Peak area, Targeted: absolute concentration (check
  "experimentalDetails" for unit)

- tissue_code:

  MoTrPAC tissue code

- tissue_name:

  Tissue name

- group_time_point:

  Intervention group (Exercise / Control) + time point

- sex:

  Animal Sex

- site_code:

  Chemical Analysis Site (CAS) short abbreviation

- group:

  Intervention group: Exercise / Control

- condition:

  Sex + group + time-point

- bioreplicate:

  Sex + group + time-point + sample_order
