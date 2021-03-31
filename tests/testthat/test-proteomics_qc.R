context("Test check_proteomics functions")

test_that("check_ratio_proteomics returns NULL when columns are missed", {
  expect_equal(check_ratio_proteomics(df_ratio = metadata_metabolites_named, isPTM =  TRUE, return_n_issues = TRUE, verbose = FALSE), 10)
})

test_that("check_rii_proteomics returns NULL when columns are missed", {
  expect_equal(check_rii_proteomics(df_rri = metadata_metabolites_named, isPTM =  TRUE, return_n_issues = TRUE), 10)
})

test_that("check_vial_metadata returns NULL when columns are missed", {
  expect_equal(check_vial_metadata_proteomics(df_vm = metadata_metabolites_named, return_n_issues = TRUE), 3)
})

