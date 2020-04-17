# R script to generate R data objects used by the library

# BIC tissue codes
bic_animal_tissue_code <- read.delim("inst/extdata/bic_animal_tissue_code.txt", stringsAsFactors = FALSE)
save(bic_animal_tissue_code, file = "data/bic_animal_tissue_code.RData")

