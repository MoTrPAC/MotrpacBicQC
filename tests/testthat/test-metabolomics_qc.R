context("Test check_metabolomics functions")

test_that("check_metadata_metabolites returns the right number of issues", {
  expect_equal(check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named", return_n_issues = TRUE, verbose = FALSE), 0)
  expect_equal(check_metadata_metabolites(df = metadata_metabolites_unnamed, name_id = "unnamed", return_n_issues = TRUE, verbose = FALSE), 0)
  expect_equal(check_metadata_metabolites(df = metadata_metabolites_unnamed, name_id = "named", return_n_issues = TRUE, verbose = FALSE), 2)
})

test_that("check_metadata_sample returns the right number of issues", {
  expect_equal(check_metadata_samples(df = metadata_sample_named, cas = "umichigan", return_n_issues = TRUE, verbose = FALSE), 3)
  expect_equal(check_metadata_samples(df = metadata_sample_unnamed, cas = "umichigan", return_n_issues = TRUE, verbose = FALSE), 3)
})

test_that("check_results returns the right number of issues", {
  expect_equal( check_results(r_m = results_named, m_m = metadata_metabolites_named, m_s = metadata_sample_named, return_n_issues = TRUE, verbose = FALSE), 0)
  expect_equal( check_results(r_m = results_unnamed, m_m = metadata_metabolites_unnamed, m_s = metadata_sample_unnamed, return_n_issues = TRUE, verbose = FALSE), 0)
  expect_equal( check_results(r_m = results_unnamed, m_m = metadata_metabolites_named, m_s = metadata_sample_unnamed, return_n_issues = TRUE, verbose = FALSE), 1)
})
