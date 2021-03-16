
# R script to generate the tissue codes and colors
library(data.table)
library(devtools)

# BIC tissue codes----
bic_animal_tissue_code <- read.delim("inst/extdata/bic_animal_tissue_code.tsv", stringsAsFactors = FALSE)
bic_animal_tissue_code = data.table(bic_animal_tissue_code)

# It adds as well tissue abbreviations and colors
bic_animal_tissue_code[,abbreviation := tissue_abbr[match(gsub(' powder','',(tolower(bic_animal_tissue_code[,bic_tissue_name]))),
                                                                names(tissue_abbr))]]
bic_animal_tissue_code[,tissue_hex_colour := tissue_cols[match(bic_animal_tissue_code[,abbreviation], names(tissue_cols))]]
# write.table(bic_animal_tissue_code, file='~/updated_bic_animal_tissue_code.tsv', sep='\t', col.names=T, row.names=F, quote=F) 

use_data(bic_animal_tissue_code, overwrite = TRUE)
