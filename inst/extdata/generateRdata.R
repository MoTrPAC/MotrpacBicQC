# R script to generate R data objects used by the library

library(dplyr)
library(devtools)

# FILES IN PACKAGE--------------------------------------------------------------

## Update `assay_codes`: add/edit

### Currently (last edit: 2025-03-13)

#' 1. Open the `assay_codes.csv` file from the inst folder of this package 
#' with your favorite editor and modify, edit, the content.
#' 
#' 2. Open the assay_codes.csv file from the inst folder of this package:
filename <- "inst/extdata/assay_codes.csv"
assay_codes <- read.csv(filename, stringsAsFactors = FALSE)

#' 3. Overwrite the data in the package with the new one:
usethis::use_data(assay_codes, overwrite = TRUE)



### Previous updates of the assay_codes (deprecated)

## update assay codes (2024-01-03)-----

### Update the "proteomics" ome by proteomics-untargeted

assay_codes <- assay_codes %>%
  mutate(omics_code = ifelse(omics_code == "proteomics", "proteomics-untargeted", omics_code))

OLINK <- data.frame(omics_code = "proteomics-targeted",
                      submission_code = "PROT_OL",
                      assay_code = "prot-ol",
                      assay_name = "Proximity extension assay-based technology for multiplexed protein analysis",
                      cas_code = "broad_rg",
                      assay_abbreviation = "OLINK",
                      assay_hex_colour = "#4DAF4A")

assay_codes <- rbind(assay_codes, OLINK)

use_data(assay_codes, overwrite = TRUE)

## update assay codes (2023-01-12)-----

IMM_GLC <- data.frame(omics_code = "metabolomics-targeted",
                        submission_code = "IMM_GLC",
                        assay_code = "metab-t-imm-glc",
                        assay_name = "Immunoassay for glucagon and related metabolic hormones",
                        cas_code = "duke",
                        assay_abbreviation = "METAB",
                        assay_hex_colour = "#E41A1C")

IMM_INS <- data.frame(omics_code = "metabolomics-targeted",
                      submission_code = "IMM_INS",
                      assay_code = "metab-t-imm-ins",
                      assay_name = "Immunoassay for insulin and related metabolic hormones",
                      cas_code = "duke",
                      assay_abbreviation = "METAB",
                      assay_hex_colour = "#E41A1C")

IMM_CRT <- data.frame(omics_code = "metabolomics-targeted",
                      submission_code = "IMM_CRT",
                      assay_code = "metab-t-imm-crt",
                      assay_name = "Immunoassay for corticosteroids",
                      cas_code = "duke",
                      assay_abbreviation = "METAB",
                      assay_hex_colour = "#E41A1C")
                     
assay_codes <- rbind(assay_codes, IMM_CRT, IMM_GLC, IMM_INS)


assay_codes$assay_name[which(assay_codes$submission_code == "CONV")] <- "Clinical chemistry assays for conventional metabolites"

use_data(assay_codes, overwrite = TRUE)


## update assay codes (code by Nicole 2023-01-12)-----
library(MotrpacBicQC)
library(MotrpacRatTraining6moData)
library(data.table)

assay_key = data.table(MotrpacBicQC::assay_codes)

# add a row for immuno 
assay_key = rbindlist(list(assay_key, 
                           data.table(omics_code="proteomics-targeted",
                                      submission_code="",
                                      assay_code="immunoassay",
                                      assay_name="Multiplexed immunoassays",
                                      cas_code="stanford")))

# add assay abbreviations
assay_key[,assay_abbreviation := ASSAY_CODE_TO_ABBREV[assay_code]]
assay_key[grepl("metab",omics_code), assay_abbreviation := "METAB"]

# add colors
assay_key[,assay_hex_colour := ASSAY_COLORS[assay_abbreviation]]

assay_codes = as.data.frame(assay_key)
use_data(assay_codes, overwrite = TRUE)

# Phenotypic data-----

# load("~/DriveStanford/motrpac/dmaqc/phenotypes.Rdata")
# phenotypes_pass1a06_short <- phenotypes[c("tissue_ecode", "vial_label", "tissue_name", "group_time_point", "sex", "site_code")]
# colnames(phenotypes_pass1a06_short)[grep("tissue_ecode", colnames(phenotypes_pass1a06_short))] <- "tissue_code"
# save(phenotypes_pass1a06_short, file = "~/github/MoTrPAC/MotrpacBicQC/data/phenotypes_pass1a06_short.RData")
# use_data(phenotypes_pass1a06_short)


# Metabolomics Workbench data dictionary

# metabolomics_data_dictionary <- read.csv("inst/extdata/motrpac-metabolomics-named-revised-20210325.csv", stringsAsFactors = FALSE)
metabolomics_data_dictionary <- read.csv("inst/extdata/motrpac-metabolomics-named-revised-20210429.csv", stringsAsFactors = FALSE)
metabolomics_data_dictionary$assay <- NULL

if(any(duplicated(metabolomics_data_dictionary$refmet_name))){
  metabolomics_data_dictionary <- metabolomics_data_dictionary[!duplicated(metabolomics_data_dictionary$refmet_name), ]  
}else{
  message("No duplicated ref_met ids")
}

use_data(metabolomics_data_dictionary, overwrite = TRUE)


# FILES NOT IN PACKAGE--------------------------------------------------------------
setwd("~/motrpac-portal-transfer-michigan/PASS1A-06/T31/IONPNEG/BATCH1_ 20190909/")

# metadata_metabolites_named, unnamed----
metadata_metabolites_named <- read.delim("PROCESSED_20200205/NAMED/metadata_metabolites_named_EX00930_PASS1A_Plasma_IP_Neg_Named_1_0_1_0_ver3.txt", stringsAsFactors = FALSE)
metadata_metabolites_unnamed_full <- read.delim("PROCESSED_20200205/UNNAMED/metadata_metabolites_unnamed_EX00930_PASS1A_Plasma_IP_Neg_Features_1_0_1_0_ver_20190927.txt", stringsAsFactors = FALSE)

# metadata_sample_named, unnamed-----
metadata_sample_named <- read.delim("PROCESSED_20200205/NAMED/metadata_sample_named_EX00930_PASS1A_Plasma_IP_Neg_Named_1_0_1_0_ver3.txt", stringsAsFactors = FALSE)
metadata_sample_named <- metadata_sample_named[c("sample_id", "sample_type", "sample_order", "raw_file")]

metadata_sample_unnamed <- read.delim("PROCESSED_20200205/UNNAMED/metadata_sample_unnamed_EX00930_PASS1A_Plasma_IP_Neg_Features_1_0_1_0_ver_20190927.txt", stringsAsFactors = FALSE)
metadata_sample_unnamed <- metadata_sample_unnamed[c("sample_id", "sample_type", "sample_order", "raw_file")]

# results_named, unamed----
results_named <- read.delim("PROCESSED_20200205/NAMED/results_metabolites_named_EX00930_PASS1A_Plasma_IP_Neg_Named_1_0_1_0_ver3.txt",
                            stringsAsFactors = FALSE, check.names = FALSE)

results_unnamed_full <- read.delim("PROCESSED_20200205/UNNAMED/results_metabolites_unnamed_EX00930_PASS1A_Plasma_IP_Neg_Features_1_0_1_0_ver_20190927.txt",
                                   stringsAsFactors = FALSE, check.names = FALSE)

# Reduce the number of unnamed metabolites to the same size than the named metabolites----
results_unnamed <- results_unnamed_full %>% dplyr::sample_n(dim(results_named)[1])

metadata_metabolites_unnamed <- metadata_metabolites_unnamed_full[which(metadata_metabolites_unnamed_full$metabolite_name %in% results_unnamed$metabolite_name),]
metadata_metabolites_named <- metadata_metabolites_named[c("metabolite_name", "refmet_name", "mz", "rt", "formula", "neutral_mass")]
metadata_metabolites_unnamed <- metadata_metabolites_unnamed[c("metabolite_name", "mz", "rt", "neutral_mass")]

# SAVE RData sets-----
# save(results_named, file = "~/github/MoTrPAC/MotrpacBicQC/data/results_named.RData")
# save(results_unnamed, file = "~/github/MoTrPAC/MotrpacBicQC/data/results_unnamed.RData")
# save(metadata_metabolites_named, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_metabolites_named.RData")
# save(metadata_metabolites_unnamed, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_metabolites_unnamed.RData")
# save(metadata_sample_named, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_sample_named.RData")
# save(metadata_sample_unnamed, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_sample_unnamed.RData")

use_data(results_named)
use_data(results_unnamed)
use_data(metadata_metabolites_named)
use_data(metadata_metabolites_unnamed)
use_data(metadata_sample_named)
use_data(metadata_sample_unnamed)


