## ----eval = FALSE-------------------------------------------------------------
#  install.packages("devtools")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("MoTrPAC/MotrpacBicQC", build_vignettes = TRUE)

## ----setup--------------------------------------------------------------------
library(MotrpacBicQC)

## ----eval = FALSE-------------------------------------------------------------
#  # Just copy and paste in the RStudio terminal
#  test <- check_ratio_proteomics(df_ratio = metadata_metabolites_named,
#                                 isPTM =  TRUE)
#  test <- check_rii_proteomics(df_rri = metadata_metabolites_named,
#                               isPTM =  TRUE)
#  test <- check_vial_metadata_proteomics(df_vm = metadata_metabolites_named)
#  

## ----eval = FALSE-------------------------------------------------------------
#  validate_proteomics(input_results_folder = "/full/path/to/RESULTS_YYYYMMDD",
#                      cas = "your_site_code")
#  
#  # in the example above...
#  validate_proteomics(input_results_folder = "/full/path/to/PASS1B-06/T55/PROT_PH/BATCH1_20200312/RESULTS_20200909",
#                      cas = "pnnl",
#                      isPTM = TRUE,
#                      return_n_issues = FALSE)
#  

## ----eval = FALSE-------------------------------------------------------------
#  # Open the ratio results file
#  
#  proteomics_ratio_results <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  
#  check_ratio_proteomics(df_ratio = proteomics_ratio_results,
#                         isPTM = TRUE,
#                         printPDF = FALSE)
#  

## ----eval = FALSE-------------------------------------------------------------
#  # Open your files
#  proteomics_ratio_results <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  
#  check_rii_proteomics(df_rri = proteomics_ratio_results, cas = "your_side_id")

## ----eval = FALSE-------------------------------------------------------------
#  # Open your files
#  vial_metadata <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  
#  check_vial_metadata_proteomics(df_vm = vial_metadata)

## ----eval=FALSE---------------------------------------------------------------
#  ?check_vial_metadata_proteomics

