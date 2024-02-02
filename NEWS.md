# MotrpacBicQC 0.9.1 (2024-02-02)

* New assays: `PROT_OX`
* Fix package conflicts
* OLINK: write release adjustments

# MotrpacBicQC 0.9.0 (2024-01-04)

* Add support for OLINK datasets (check `olink_qc` vignette to find out more)
* Adjust function to download data from GCP (`dl_read_gcp`): 
it automatically detects the operating system (arguments `ignore_std_err` and
`ignore_std_out` deprecated)
* Multiple fixes and enhancements


# MotrpacBicQC 0.8.9 (2023-07-07)

* Fix bug preventing the processing of BICRESULTS folders (proteomics)
* Make clear that the metadata_phase.txt file is required
* Other enhancements

# MotrpacBicQC 0.8.8 (2023-06-15)

* Add 24 hours time support for the `acquisition_date` (`MM/DD/YYYY HH:MM:SS`)

# MotrpacBicQC 0.8.7 (2023-05-22)

* Add QC for the new required batching variables
* Replace deprecated `ggplot` function
* Fix issues with `dl_read_gcp`
* Other adjustments

# MotrpacBicQC 0.8.6 (2023-04-21)

* Minor adjustments

# MotrpacBicQC 0.8.5 (2023-03-18)

* New tissue codes, abbreviations, and colors available for lateral Gastrocnemius and vena cava

# MotrpacBicQC 0.8.4 (2023-02-20)

* Add new `dl_read_gcp`
* Replace dplyr `summarise` function (deprecated) by `reframe`
* Fix bug affecting `proteomics_plots`

# MotrpacBicQC 0.8.3 (2023-02-08)

* Fix typo

# MotrpacBicQC 0.8.1 and 0.8.2 (2023-02-08)

* Fix bug affecting IMM assays

# MotrpacBicQC 0.8.0 (2023-02-05)

* New metabolomics targeted assays (IMM_GLC, IMM_INS, IMM_CTR)
* Remove exception of non-unique raw files allowed for CONV assays, and added to IMM assays
* Improve metabolomics documentation


# MotrpacBicQC 0.7.9 (2023-01-12)

* Add exception: unique raw files are now not required for metabolomics CONV assay
* Update `assay_codes`: immunoassay/IMMUNO added. The table now also includes assay hex colours and assay abbreviation

# MotrpacBicQC 0.7.8 (2022-12-13)

* Bug fix

# MotrpacBicQC 0.7.7 (2022-12-08)

* Better handling of large proteomics datasets:
   * Proteomics RII plots are skipped if the dataset is too large
   * Larger pdf size for proteomics ratio plots
* Several improvements and enhancements

# MotrpacBicQC 0.7.6 (2022-10-27)

* Improve DMAQC validation

# MotrpacBicQC 0.7.5 (2022-10-21)

* Update `write_proteomics` according to latest updates on data/file structure
* Fix bugs affecting `metabolomics_qc` and checks on `file_manifest`

# MotrpacBicQC 0.7.4 (2022-10-09)

* Metabolomics plots: check if enough compounds to generate plots
* Update `assay_codes`: conventional assays code added (CONV)

# MotrpacBicQC 0.7.3 (2022-10-02)

* Adjustments to generate data releases (deal with pass1a/1c)
* Version number will be added to upcoming releases

# MotrpacBicQC 0.7.2 (2022-09-22)

* Update a package's dependency

# MotrpacBicQC 0.7.1 (2022-09-17)

Updates affecting the proteomics validation:

* `validate_proteomics`: rename argument `run_by_bic` to `check_only_results`. Default is still `FALSE` (it does not affect CAS)
* Adjust size of pdf output depending on the number of samples
* Support `BICRESULTS_YYYYMMDD` folder validation (similar to the currently supported `RESULTS_YYYYMMDD` and `PROCESSED_YYYYMMDD` folders). This folder is the output of the proteomics pipeline run by the BIC

# MotrpacBicQC 0.7.0 (2022-08-25)

* Update motrpac color abbreviations

# MotrpacBicQC 0.6.9 (2022-07-09)

* Metabolomics: new density plots
* Code optimizations
* Bug fixes affecting file manifest checks
   
   
# MotrpacBicQC 0.6.8 (2022-04-11)

* Refactor the DMAQC validation. A new file will be required when:
   + Two phases are combined in the same batch (e.g., `PASS1A-06|PASS1C-06`)
   + The phase content is different from the input folder name (e.g., `PASS1C-06` might 
   be submitted but the input folder name is `PASS1A-06`)

# MotrpacBicQC 0.6.7 (2022-03-11)

* DMAQC validation: print out missing vial labels
* Metabolomics QC: add mz/rt density plots

# MotrpacBicQC 0.6.6 (2022-03-03)

* Improve reporting and handling of required files
* Fix minor bugs 

# MotrpacBicQC 0.6.5 (2022-02-28)

* Metabolomics: fix manifest checks.

# MotrpacBicQC 0.6.4 (2022-02-27)

* Metabolomics: adjust metabolomic plots to deal with a large number of samples
* Metabolomics new plot: sum of intensity/concentration
* Metabolomics: detects negative values
* Metabolomics: update vignette

# MotrpacBicQC 0.6.3 (2022-02-25)

* Proteomics: support for tmt16
* Metabolomics: Improve verbosity for wrong tissue code

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
