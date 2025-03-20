## ----eval = FALSE-------------------------------------------------------------
# install.packages("devtools")

## ----eval = FALSE-------------------------------------------------------------
# library(devtools)
# devtools::install_github("MoTrPAC/MotrpacBicQC", build_vignettes = TRUE)

## ----setup--------------------------------------------------------------------
library(MotrpacBicQC)

## ----eval = FALSE-------------------------------------------------------------
# # Just copy and paste in the RStudio terminal.
# check_metadata_samples_lab(df = metadata_metabolites_named)
# check_metadata_analyte(df = metadata_metabolites_named)
# check_results_assays(df = results_named, assay_type = "lab")

## ----eval = FALSE-------------------------------------------------------------
# n_issues <- validate_lab(input_results_folder = "/full/path/to/HUMAN/T02/LAB_CK/BATCH1_20221102/PROCESSED_20221102/",
#                          cas = "duke",
#                          return_n_issues = TRUE,
#                          verbose = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# ?validate_lab

