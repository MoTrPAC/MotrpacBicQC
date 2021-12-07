# MotrpacBicQC 0.6.2 (2021-12-07)

* Fix bug affecting the validation of refmet_names
* Update data objects (immunoassay added)

# MotrpacBicQC 0.6.1 (2021-09-06)

* Refactor the validation of refmet_name. It now checks on at the time 
using the RefMet API. It also validates multipeak isoforms
* The function `get_and_validate_mdd()` donwload the entire RefMet database (warning, >15MB)
* The `metabolomics_data_dictionary` data object will be soon deprecated. 

# MotrpacBicQC 0.6.0 (2021-09-06)

* Support DMAQC validation of human submissions

# MotrpacBicQC 0.5.9 (2021-08-10)

* Adjustments for PASS1C-06
* Bug fixes affecting HUMAN phase processing

# MotrpacBicQC 0.5.8 (2021-08-06)

* New Metabolomics QC plots: number and proportion of named vs unnamed features identified

# MotrpacBicQC 0.5.7 (2021-08-05)

* New proteomics QC plots for protein coverage

# MotrpacBicQC 0.5.6 (2021-07-22)

* Bug fix

# MotrpacBicQC 0.5.5 (2021-07-21)

* Update: add human tissue codes

# MotrpacBicQC 0.5.4 (2021-07-16)

* Improved version for checking the file manifest from metabolomics submissions
* Enable DMAQC validation for submissions combining multiple phases (e.g. PASS1A-06 + PASS1C-06)

# MotrpacBicQC 0.5.3 (2021-06-25)

* Metabolomics QC: new metabolomics QC plots, including number of IDs per sample, intensity distribution, and percentage of NA values
* Markdown: Replace `prettydoc` by `rmdformats`
* New assay code: CONV (Targeted Conventional metabolites or clinical analytes, provided by Duke)

# MotrpacBicQC 0.5.2 (2021-06-23)

* New Phase: HUMAN (name of the new project folder for the human studies)
* Proteomics QC: new proteomics QC plot, number of unique ids per sample
* Proteomics QC: improved QC plots

# MotrpacBicQC 0.5.1 (2021-05-03)

* New metabolomics `sample_type`: `QC-ReCAS`, Global reference biological material prepared at

# MotrpacBicQC 0.5.0 (2021-04-29)

* Update data dictionary (GTech's kegg revision)

# MotrpacBicQC 0.4.9 (2021-04-27)

* Bug fixes
* Update readme

# MotrpacBicQC 0.4.8 (2021-04-20)

* Proteomics QC: enable option to check data processed at the BIC
* Bug fixes (color code)


# MotrpacBicQC 0.4.7 (2021-04-02)

* Proteomics QC: update warnings affecting `gene_symbol` and `entrez_id` when missing ids
* Fix issue affecting Windows machines (Pierre J-B)
* Several Color code fixes * assay (Nicole G)
* Bug fixes

# MotrpacBicQC 0.4.6 (2021-03-26)

* New QC checks for the new proteomics requirements

# MotrpacBicQC 0.4.5 (2021-03-25)

* Update metabolomics data dictionary, including:
   * Broad Metabolomics revision (39 new kegg ids * minor corrections)
   * targeted refmet_name 
* Fix text in warning in case of missing manifest file
* Add and improve tests

# MotrpacBicQC 0.4.4 (2021-03-22)

* Fix and improve check (new required)
 manifest file for untargeted metabolomics datasets

# MotrpacBicQC 0.4.3 (2021-03-16)

* Update `tissue_cols` data object (Nicole Gay)


# MotrpacBicQC 0.4.2 (2021-03-10)

* Fix bug in `validate_metabolomics` affecting unnamed sample check

# MotrpacBicQC 0.4.1 (2021-03-08)

* Improved QC for the required manifest file (proteomics and metabolomics)


# MotrpacBicQC 0.4.0 (2021-03-04)

* Metabolomics data dictionary available as a data object: `metabolomics_data_dictionary`
* The function `get_and_validate_mdd()`
 still works, but does not pull the data from Metabolomics Workbench. Just return metabolomics_data_dictionary
* Bug fix: the manifest issue count now properly displays the number of issues detected.

# MotrpacBicQC 0.3.9 (2021-02-10)

* Bug fix: restore previous version of `bic_animal_tissue_code`

# MotrpacBicQC 0.3.8 (2021-02-03)

* Add new assays

# MotrpacBicQC 0.3.5 (2020-10-20)

* Bug fix on DMAQC to deal with missing data in DMAQC table
* Proteomics: Address 130C missing channel issue affecting the Broad
* Proteomics: add an additional QC check point to validate that all values in vial_label column are unique
* Proteomics: `write_proteomics` updated (only required columns selected)
* Fix bugs and typos

# MotrpacBicQC 0.3.4 (2020-10-15)

* Add a check point when expected files are not available in the manifest file

# MotrpacBicQC 0.3.3 (2020-10-14)

* Fix bugs affecting manifest files
* Adjustments when files don't meet requirements

# MotrpacBicQC 0.3.2 (2020-10-06)

* Check new required manifest file in both proteomics and metabolomics submissions
* Raw file manifest now optional

# MotrpacBicQC 0.3.1 (2020-10-04)

* Proteomics write proteomics
* Proteomics load metabolomics

# MotrpacBicQC 0.3.0 (2020-09-21)

* Proteomics QC support

# MotrpacBicQC 0.2.2 (2020-08-09)

* New assays

# MotrpacBicQC 0.2.1 (2020-06-03)

* Bug fixes
