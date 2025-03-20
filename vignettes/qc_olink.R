## ----eval = FALSE-------------------------------------------------------------
#  install.packages("devtools")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("MoTrPAC/MotrpacBicQC", build_vignettes = TRUE)

## ----setup--------------------------------------------------------------------
library(MotrpacBicQC)

## ----eval = FALSE-------------------------------------------------------------
#  # Just copy and paste in the RStudio terminal.
#  check_metadata_samples_olink(df = metadata_metabolites_named)
#  check_metadata_proteins(df = metadata_metabolites_named)
#  check_results_assays(df = results_named, assay_type = "olink")

## ----eval = FALSE-------------------------------------------------------------
#  n_issues <- validate_olink(input_results_folder = "/full/path/to/HUMAN/T02/PROT_OL/BATCH1_20210825/RESULTS_20221102",
#                             cas = "broad_rg",
#                             return_n_issues = TRUE,
#                             verbose = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  ?validate_olink

