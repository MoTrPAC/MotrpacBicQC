## ----eval = FALSE-------------------------------------------------------------
#  install.packages("devtools")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("MoTrPAC/MotrpacBicQC", build_vignettes = TRUE)

## ----setup--------------------------------------------------------------------
library(MotrpacBicQC)

## ----eval = FALSE-------------------------------------------------------------
#  # Just copy and paste in the RStudio terminal
#  
#  check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named")
#  check_metadata_samples(df = metadata_sample_named, cas = "umichigan")
#  check_results(r_m = results_named, m_s = metadata_sample_named, m_m = metadata_metabolites_named)

## ----eval = TRUE--------------------------------------------------------------
check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named")
check_metadata_samples(df = metadata_sample_named, cas = "umichigan")
check_results(r_m = results_named, m_s = metadata_sample_named, m_m = metadata_metabolites_named)

## ----eval = FALSE-------------------------------------------------------------
#  validate_metabolomics(input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD",
#                        cas = "your_site_code")

## ----eval = FALSE-------------------------------------------------------------
#  # Open the metadata_metabolites file(s)
#  
#  metadata_metabolites_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  metadata_metabolites_unnamed <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  
#  check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named")
#  check_metadata_metabolites(df = metadata_metabolites_unnamed, name_id = "unnamed")
#  

## ----eval = FALSE-------------------------------------------------------------
#  # Open your files
#  metadata_sample_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  metadata_sample_unnamed <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  
#  check_metadata_samples(df = metadata_sample_named, cas = "your_side_id")
#  check_metadata_samples(df = metadata_sample_unnamed, cas = "your_side_id")

## ----eval = FALSE-------------------------------------------------------------
#  # Open your files
#  metadata_metabolites_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  metadata_sample_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  results_named <- read.delim(file = "/path/to/your/file", stringsAsFactors = FALSE)
#  
#  check_results(r_m = results_named,
#                m_s = metadata_sample_named,
#                m_m = metadata_metabolites_named)

## ---- eval=FALSE--------------------------------------------------------------
#  t31_ionpneg <- combine_metabolomics_batch(input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD/",
#                                            cas = "umichigan")

## ---- eval=FALSE--------------------------------------------------------------
#  plasma.untargeted.merged <-
#    merge_all_metabolomics(m_m_n = metadata_metabolites_named,
#                           m_m_u = metadata_metabolites_unnamed,
#                           m_s_n = metadata_sample_named,
#                           r_n = results_named,
#                           r_u = results_unnamed,
#                           phase = "PASS1A-06")

## ----eval=FALSE---------------------------------------------------------------
#  ?merge_all_metabolomics

