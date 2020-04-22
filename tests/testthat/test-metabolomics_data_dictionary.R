context("Test metabolomics data dictionary")


test_that("REST service from RefMet works and return right number of columns and values", {
  expected_mdd_colnames <- c("refmet_name", "super_class", "main_class", "sub_class", "formula", "exactmass", "inchi_key", "pubchem_cid", "standard")
  mdd <- get_and_validate_mdd()
  expect_equal( setequal(expected_mdd_colnames, colnames(mdd)), TRUE)
  expect_gt(dim(mdd)[1], 2000)
})
