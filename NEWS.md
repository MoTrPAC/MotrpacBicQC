
# MotrpacBicQC 1.3.0 (2025-12-06)

## New Features

* **Human Study Support**: Enhanced handling of HUMAN phase formats, including main and pre-COVID tranches, for more accurate validation of human study submissions

## Bug Fixes

* **DMAQC Validation**: Fixed a critical bug where DMAQC validation failures were not being counted in the total issues. The `check_viallabel_dmaqc` function returns a string status (`"OK"`, `"FAIL"`), which was incorrectly checked with `is.numeric()`. Now properly increments the issue count when validation fails
* **Proteomics**: Fixed missing sample labels when reference (Ref) channels are absent in proteomics datasets
* **Metabolomics**: Fixed incorrect variable assignment in `write_metabolomics_releases` where cleaned sample metadata was incorrectly assigned to `m_s_u` instead of `m_s_n`
* **LAB QC Plots**: Fixed plot layout issues in `plot_basic_lab_qc` that caused plots to print to console instead of saving to PDF. Refactored to use `gridExtra::grid.arrange()` properly

## Improvements

* **QC Date Format**: Standardized QC date format to `YYYYMMDD` across all validation modules (`validate_lab`, `validate_metabolomics`, `validate_olink`, `validate_proteomics`), ensuring consistency in output filenames and reports
* **Phase Parsing**: Improved phase parsing in `check_viallabel_dmaqc` to handle complex HUMAN phase formats (e.g., `HUMAN-MAIN-TR04`) correctly
* **Release Writer**: Added handling for `human-main` phase in release folder logic for proper directory structure generation
* **PDF Output**: Improved plot output flow—suppresses console noise when saving to PDF and logs saved file locations

## Internal Changes
* Updated `.Rbuildignore` to use explicit regex patterns for more reliable build process
* Various code refactoring for better maintainability

---

# MotrpacBicQC 1.2.0 (2025-07-30)

## New Features

* **Assay Codes**: Added new assay code for Whole Genome Sequencing (WGS) to `inst/extdata/assay_codes.csv`

## Improvements

* **Assay Codes Updates**: Populated `submission_code` values for several existing transcriptomics and epigenomics entries; corrected typo for 'Phosphoproteomics'
* **File Location Enforcement**: Modified `set_phase` function to enforce that `metadata_phase.txt` must be located directly within the batch folder (changed `recursive` parameter from `TRUE` to `FALSE`)

---
  
# MotrpacBicQC 1.1.0 (2025-05-01)

## Improvements

* Updated the data object `assay_codes`:
  + Updated `omics_text` values to omics technologies
  + Added `ome_text` for the "omes": the complete set of a given biological entity (all genes, all proteins, etc.)
  + Added `assay` as a copy of `assay_code`

---
  
# MotrpacBicQC 1.0.0 (2025-04-23)

## New Features

* Updated the data object `assay_codes` with 2 new variables:
  + `assay_short_text`: Abbreviated name specific to the assay for use in graphs and tables
  + `ome_text`: Omic class measured by the given assay

## Improvements

* Updated style of the vignettes

---

# MotrpacBicQC 0.9.9 (2025-04-22)

* Fixed minor issue with `dl_read_gcp`

---

# MotrpacBicQC 0.9.8 (2025-04-10)

* Updated `dl_read_gcp`: now supports `gcloud` in addition to `gsutil`

---

# MotrpacBicQC 0.9.7 (2025-03-16)

## New Features

* **Clinical Chemistry Support**: Added QC support for clinical chemistry assays: glucagon, insulin, cortisol, and creatine kinase
* Conventional metabolites (previously `metab-t-conv`) now expected as a new assay within this category (`lab-conv`)

## Bug Fixes

* Numerous bug fixes and enhancements

---

# MotrpacBicQC 0.9.6 (2024-09-23)

* Download and read file from GCP function can create recursive folders (@christopherjin)
* Adjustments in metabolomics metadata sample files QC to enable processing of old submissions (before batch related variables were required)

---

# MotrpacBicQC 0.9.5 (2024-05-22)

* **Proteomics**: Added QC support for TMT-18

---
  
# MotrpacBicQC 0.9.4 (2024-05-16)

* Enhanced and improved `dl_read_gcp`: 
  + Check if `gsutil` path is correct and report back to the user if it is not
  + Handle spaces in folder names (although not recommended)
  + Improved error source detection
  + Improved verbosity and feedback to the user

---
  
# MotrpacBicQC 0.9.3 (2024-03-25)

* **Critical Update**: Fixed `validate_refmetname` to ensure checking the refmet standardized name; updated refmet tests
* Updated `get_and_validate_mdd()`:
  + Updated REST service URL
  + Updated documentation
  + Removed dependency on data.table
* Enhanced: only one `metadata_phase` file allowed
* Enhanced `dl_read_gcp`: replaced data.table by read_delim
* Enhanced `open_file`: accepts only tab-delimited files

---

# MotrpacBicQC 0.9.2 (2024-03-04)

* **Critical Update**: Resolved an issue where the validation of refmet names was compromised due to updates to the Metabolomics Workbench REST service. This version introduces adjustments to ensure accurate validation of refmet names.

---

# MotrpacBicQC 0.9.1 (2024-02-02)

* New assay: `PROT_OX`
* Fixed package conflicts
* OLINK: write release adjustments

---

# MotrpacBicQC 0.9.0 (2024-01-04)

## New Features

* Added support for OLINK datasets (check `olink_qc` vignette to find out more)

## Improvements

* Adjusted function to download data from GCP (`dl_read_gcp`): automatically detects the operating system (arguments `ignore_std_err` and `ignore_std_out` deprecated)
* Multiple fixes and enhancements

---

# MotrpacBicQC 0.8.9 (2023-07-07)

* Fixed bug preventing the processing of BICRESULTS folders (proteomics)
* Made clear that the `metadata_phase.txt` file is required
* Other enhancements

---

# MotrpacBicQC 0.8.8 (2023-06-15)

* Added 24-hour time support for the `acquisition_date` (`MM/DD/YYYY HH:MM:SS`)

---

# MotrpacBicQC 0.8.7 (2023-05-22)

* Added QC for the new required batching variables
* Replaced deprecated `ggplot` function
* Fixed issues with `dl_read_gcp`
* Other adjustments

---

# MotrpacBicQC 0.8.6 (2023-04-21)

* Minor adjustments

---

# MotrpacBicQC 0.8.5 (2023-03-18)

* New tissue codes, abbreviations, and colors available for lateral Gastrocnemius and vena cava

---

# MotrpacBicQC 0.8.4 (2023-02-20)

* Added new `dl_read_gcp`
* Replaced dplyr `summarise` function (deprecated) by `reframe`
* Fixed bug affecting `proteomics_plots`

---

# MotrpacBicQC 0.8.3 (2023-02-08)

* Fixed typo

---

# MotrpacBicQC 0.8.1 and 0.8.2 (2023-02-08)

* Fixed bug affecting IMM assays

---

# MotrpacBicQC 0.8.0 (2023-02-05)

## New Features

* New metabolomics targeted assays (IMM_GLC, IMM_INS, IMM_CTR)

## Bug Fixes

* Removed exception of non-unique raw files allowed for CONV assays, and added to IMM assays

## Improvements

* Improved metabolomics documentation

---

# MotrpacBicQC 0.7.9 (2023-01-12)

* Added exception: unique raw files are now not required for metabolomics CONV assay
* Updated `assay_codes`: immunoassay/IMMUNO added. The table now also includes assay hex colours and assay abbreviation

---

# MotrpacBicQC 0.7.8 (2022-12-13)

* Bug fix

---

# MotrpacBicQC 0.7.7 (2022-12-08)

* Better handling of large proteomics datasets:
   * Proteomics RII plots are skipped if the dataset is too large
   * Larger PDF size for proteomics ratio plots
* Several improvements and enhancements

---

# MotrpacBicQC 0.7.6 (2022-10-27)

* Improved DMAQC validation

---

# MotrpacBicQC 0.7.5 (2022-10-21)

* Updated `write_proteomics` according to latest updates on data/file structure
* Fixed bugs affecting `metabolomics_qc` and checks on `file_manifest`

---

# MotrpacBicQC 0.7.4 (2022-10-09)

* Metabolomics plots: check if enough compounds to generate plots
* Updated `assay_codes`: conventional assays code added (CONV)

---

# MotrpacBicQC 0.7.3 (2022-10-02)

* Adjustments to generate data releases (deal with pass1a/1c)
* Version number will be added to upcoming releases

---

# MotrpacBicQC 0.7.2 (2022-09-22)

* Updated a package's dependency

---

# MotrpacBicQC 0.7.1 (2022-09-17)

Updates affecting the proteomics validation:

* `validate_proteomics`: renamed argument `run_by_bic` to `check_only_results`. Default is still `FALSE` (does not affect CAS)
* Adjusted size of PDF output depending on the number of samples
* Support `BICRESULTS_YYYYMMDD` folder validation (similar to the currently supported `RESULTS_YYYYMMDD` and `PROCESSED_YYYYMMDD` folders). This folder is the output of the proteomics pipeline run by the BIC

---

# MotrpacBicQC 0.7.0 (2022-08-25)

* Updated MoTrPAC color abbreviations

---

# MotrpacBicQC 0.6.9 (2022-07-09)

* Metabolomics: new density plots
* Code optimizations
* Bug fixes affecting file manifest checks

---
   
# MotrpacBicQC 0.6.8 (2022-04-11)

* Refactored the DMAQC validation. A new file will be required when:
   + Two phases are combined in the same batch (e.g., `PASS1A-06|PASS1C-06`)
   + The phase content is different from the input folder name (e.g., `PASS1C-06` might be submitted but the input folder name is `PASS1A-06`)

---

# MotrpacBicQC 0.6.7 (2022-03-11)

* DMAQC validation: print out missing vial labels
* Metabolomics QC: added mz/rt density plots

---

# MotrpacBicQC 0.6.6 (2022-03-03)

* Improved reporting and handling of required files
* Fixed minor bugs 

---

# MotrpacBicQC 0.6.5 (2022-02-28)

* Metabolomics: fixed manifest checks

---

# MotrpacBicQC 0.6.4 (2022-02-27)

* Metabolomics: adjusted metabolomic plots to deal with a large number of samples
* Metabolomics new plot: sum of intensity/concentration
* Metabolomics: detects negative values
* Metabolomics: updated vignette

---

# MotrpacBicQC 0.6.3 (2022-02-25)

* Proteomics: support for TMT-16
* Metabolomics: improved verbosity for wrong tissue code

---

# MotrpacBicQC 0.6.2 (2021-12-07)

* Fixed bug affecting the validation of refmet_names
* Updated data objects (immunoassay added)

---

# MotrpacBicQC 0.6.1 (2021-09-06)

* Refactored the validation of refmet_name. It now checks one at a time using the RefMet API. It also validates multipeak isoforms
* The function `get_and_validate_mdd()` downloads the entire RefMet database (warning: >15MB)
* The `metabolomics_data_dictionary` data object will be deprecated soon

---

# MotrpacBicQC 0.6.0 (2021-09-06)

* Support DMAQC validation of human submissions

---

# MotrpacBicQC 0.5.9 (2021-08-10)

* Adjustments for PASS1C-06
* Bug fixes affecting HUMAN phase processing

---

# MotrpacBicQC 0.5.8 (2021-08-06)

* New Metabolomics QC plots: number and proportion of named vs unnamed features identified

---

# MotrpacBicQC 0.5.7 (2021-08-05)

* New proteomics QC plots for protein coverage

---

# MotrpacBicQC 0.5.6 (2021-07-22)

* Bug fix

---

# MotrpacBicQC 0.5.5 (2021-07-21)

* Added human tissue codes

---

# MotrpacBicQC 0.5.4 (2021-07-16)

* Improved version for checking the file manifest from metabolomics submissions
* Enabled DMAQC validation for submissions combining multiple phases (e.g. PASS1A-06 + PASS1C-06)

---

# MotrpacBicQC 0.5.3 (2021-06-25)

* Metabolomics QC: new metabolomics QC plots, including number of IDs per sample, intensity distribution, and percentage of NA values
* Markdown: replaced `prettydoc` by `rmdformats`
* New assay code: CONV (Targeted Conventional metabolites or clinical analytes, provided by Duke)

---

# MotrpacBicQC 0.5.2 (2021-06-23)

* New Phase: HUMAN (name of the new project folder for the human studies)
* Proteomics QC: new proteomics QC plot, number of unique IDs per sample
* Proteomics QC: improved QC plots

---

# MotrpacBicQC 0.5.1 (2021-05-03)

* New metabolomics `sample_type`: `QC-ReCAS`, Global reference biological material prepared at CAS

---

# MotrpacBicQC 0.5.0 (2021-04-29)

* Updated data dictionary (GTech's KEGG revision)

---

# MotrpacBicQC 0.4.9 (2021-04-27)

* Bug fixes
* Updated README

---

# MotrpacBicQC 0.4.8 (2021-04-20)

* Proteomics QC: enabled option to check data processed at the BIC
* Bug fixes (color code)

---

# MotrpacBicQC 0.4.7 (2021-04-02)

* Proteomics QC: updated warnings affecting `gene_symbol` and `entrez_id` when missing IDs
* Fixed issue affecting Windows machines (Pierre J-B)
* Several color code fixes for assay (Nicole G)
* Bug fixes

---

# MotrpacBicQC 0.4.6 (2021-03-26)

* New QC checks for the new proteomics requirements

---

# MotrpacBicQC 0.4.5 (2021-03-25)

* Updated metabolomics data dictionary, including:
   * Broad Metabolomics revision (39 new KEGG IDs + minor corrections)
   * Targeted refmet_name 
* Fixed text in warning in case of missing manifest file
* Added and improved tests

---

# MotrpacBicQC 0.4.4 (2021-03-22)

* Fixed and improved check for new required manifest file for untargeted metabolomics datasets

---

# MotrpacBicQC 0.4.3 (2021-03-16)

* Updated `tissue_cols` data object (Nicole Gay)

---

# MotrpacBicQC 0.4.2 (2021-03-10)

* Fixed bug in `validate_metabolomics` affecting unnamed sample check

---

# MotrpacBicQC 0.4.1 (2021-03-08)

* Improved QC for the required manifest file (proteomics and metabolomics)

---

# MotrpacBicQC 0.4.0 (2021-03-04)

* Metabolomics data dictionary available as a data object: `metabolomics_data_dictionary`
* The function `get_and_validate_mdd()` still works, but does not pull the data from Metabolomics Workbench. Just returns `metabolomics_data_dictionary`
* Bug fix: the manifest issue count now properly displays the number of issues detected

---

# MotrpacBicQC 0.3.9 (2021-02-10)

* Bug fix: restored previous version of `bic_animal_tissue_code`

---

# MotrpacBicQC 0.3.8 (2021-02-03)

* Added new assays

---

# MotrpacBicQC 0.3.5 (2020-10-20)

* Bug fix on DMAQC to deal with missing data in DMAQC table
* Proteomics: addressed 130C missing channel issue affecting the Broad
* Proteomics: added additional QC check point to validate that all values in vial_label column are unique
* Proteomics: `write_proteomics` updated (only required columns selected)
* Fixed bugs and typos

---

# MotrpacBicQC 0.3.4 (2020-10-15)

* Added check point when expected files are not available in the manifest file

---

# MotrpacBicQC 0.3.3 (2020-10-14)

* Fixed bugs affecting manifest files
* Adjustments when files don't meet requirements

---

# MotrpacBicQC 0.3.2 (2020-10-06)

* Check new required manifest file in both proteomics and metabolomics submissions
* Raw file manifest now optional

---

# MotrpacBicQC 0.3.1 (2020-10-04)

* Proteomics write proteomics
* Proteomics load metabolomics

---

# MotrpacBicQC 0.3.0 (2020-09-21)

* Proteomics QC support

---

# MotrpacBicQC 0.2.2 (2020-08-09)

* New assays

---

# MotrpacBicQC 0.2.1 (2020-06-03)

* Bug fixes
