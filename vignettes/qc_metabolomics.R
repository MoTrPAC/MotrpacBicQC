## ----eval = FALSE-------------------------------------------------------------
#  install.packages("devtools")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("MoTrPAC/MotrpacBicQC", build_vignettes = FALSE)

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
#  validate_metabolomics(input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD",
#                        cas = "your_site_code",
#                        f_proof = TRUE,
#                        out_qc_folder = "/path/to/the/folder/to/save/plots/",
#                        printPDF = TRUE)

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

## ----eval=TRUE----------------------------------------------------------------
?validate_metabolomics

