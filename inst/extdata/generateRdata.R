# R script to generate R data objects used by the library

library(dplyr)

# FILES IN PACKAGE--------------------------------------------------------------
# BIC tissue codes----
bic_animal_tissue_code <- read.delim("inst/extdata/bic_animal_tissue_code.txt", stringsAsFactors = FALSE)
save(bic_animal_tissue_code, file = "data/bic_animal_tissue_code.RData")

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
save(results_named, file = "~/github/MoTrPAC/MotrpacBicQC/data/results_named.RData")
save(results_unnamed, file = "~/github/MoTrPAC/MotrpacBicQC/data/results_unnamed.RData")
save(metadata_metabolites_named, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_metabolites_named.RData")
save(metadata_metabolites_unnamed, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_metabolites_unnamed.RData")
save(metadata_sample_named, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_sample_named.RData")
save(metadata_sample_unnamed, file = "~/github/MoTrPAC/MotrpacBicQC/data/metadata_sample_unnamed.RData")

