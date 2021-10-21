context("Test metabolomics data dictionary")


test_that("Metabolomics data dictionary", {
  # expected_mdd_colnames <- c("CURRENT_REFMET_NAME", "refmet_name", "metabolite_name", "is_standard", "super_class", "main_class", "sub_class", "formula", "exactmass", "pubchem_cid", "kegg_id", "inchi_key", "lm_id", "hmdb_id", "chebi_id")
  expected_mdd_colnames <- c("refmet_name")
  mdd <- get_and_validate_mdd()
  expect_equal( setequal(expected_mdd_colnames, colnames(mdd)), TRUE)
  expect_gt(dim(mdd)[1], 2000)
  expect_false(any(duplicated(mdd$refmet_name)))
})

# test_that("REST service from RefMet works and return right number of columns and values", {
#   expected_mdd_colnames <- c("refmet_name", "super_class", "main_class", "sub_class", "formula", "exactmass", "inchi_key", "pubchem_cid", "kegg_id", "standard")
#   mdd <- get_and_validate_mdd()
#   expect_equal( setequal(expected_mdd_colnames, colnames(mdd)), TRUE)
#   expect_gt(dim(mdd)[1], 2000)
# })
