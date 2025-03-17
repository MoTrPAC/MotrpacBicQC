context("Test check_metabolomics functions")

test_that("check_metadata_metabolites returns the right number of issues", {
  expect_equal(check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named", return_n_issues = TRUE, verbose = TRUE), 1)
  expect_equal(check_metadata_metabolites(df = metadata_metabolites_unnamed, name_id = "unnamed", return_n_issues = TRUE, verbose = FALSE), 0)
  expect_equal(check_metadata_metabolites(df = metadata_metabolites_unnamed, name_id = "named", return_n_issues = TRUE, verbose = FALSE), 2)
})

test_that("check_metadata_sample returns the right number of issues", {
  expect_equal(check_metadata_samples(df = metadata_sample_named, cas = "umichigan", return_n_issues = TRUE, verbose = FALSE), 2)
  expect_equal(check_metadata_samples(df = metadata_sample_unnamed, cas = "umichigan", return_n_issues = TRUE, verbose = FALSE), 2)
})

test_that("check_results returns the right number of issues", {
  expect_equal( check_results(r_m = results_named, m_m = metadata_metabolites_named, m_s = metadata_sample_named, return_n_issues = TRUE, verbose = FALSE), 0)
  expect_equal( check_results(r_m = results_unnamed, m_m = metadata_metabolites_unnamed, m_s = metadata_sample_unnamed, return_n_issues = TRUE, verbose = FALSE), 0)
  expect_equal( check_results(r_m = results_unnamed, m_m = metadata_metabolites_named, m_s = metadata_sample_unnamed, return_n_issues = TRUE, verbose = FALSE), 1)
})

test_that("Check that validate_dates_times function works correctly", {
  # Notice that there are many WRONG dates here
  df <- data.frame(id = 1:6,
                   datetime = c("12/31/2023 11:59 PM", "01/01/2024 00:00 AM", 
                                "02/29/2024 12:00 PM", "13/01/2024 01:00 PM", 
                                "02/28/2024 24:00 PM", "02/30/2024 12:00 PM"))
  
  expect_equal(MotrpacBicQC::validate_dates_times(df, "datetime", verbose = FALSE), 1)
  
  # Now everything fixed:
  df$datetime <- c("12/31/2023 11:59 PM", "01/01/2024 00:00 AM", 
                   "02/28/2024 12:00 PM", "12/01/2024 01:00 PM", 
                   "02/25/2024 23:00 PM", "02/27/2024 12:00 PM")
  
  expect_equal(MotrpacBicQC::validate_dates_times(df, "datetime", verbose = FALSE), 0)
  
  # Test that the function returns 0 when all date-times are in the correct format
  df$datetime <- c("01/31/2024 11:59:59 PM", "02/01/2024 12:00:00 AM", 
                   "02/29/2024 12:00:00 PM", "03/01/2024 01:00:00 PM", 
                   "02/28/2024 12:00:00 PM", "")
  
  expect_equal(MotrpacBicQC::validate_dates_times(df, "datetime", verbose = FALSE), 2)
  
  df$datetime <- c("01/31/2024 11:59:59 PM", "02/01/2024 12:00:00 AM", 
                   "02/29/2024 12:00:00 PM", "03/01/2024 01:00:00 PM", 
                   NA, "")
  
  expect_equal(MotrpacBicQC::validate_dates_times(df, "datetime", verbose = FALSE), 3)
  
  expect_error(MotrpacBicQC::validate_dates_times(df, "nonexistent_column", verbose = FALSE))
})

test_that("Check that validate_lc_column_id function works correctly", {
  # Data frame for testing
  df <- data.frame(id = 1:5,
                   lc_column_id = c("ABC", " DEF", "GHI", "JKL", "MNO"))
  
  expect_equal(MotrpacBicQC::validate_lc_column_id(df, "lc_column_id", verbose = FALSE), 1)
  
  df$lc_column_id <- c("ABC", "DEF", "GHI", "JKL", NA)
  
  expect_equal(MotrpacBicQC::validate_lc_column_id(df, column_name = "lc_column_id", verbose = FALSE), 1)
  
  expect_error(MotrpacBicQC::validate_lc_column_id(df = df, column_name = "nonexistent_column", verbose = FALSE))
})

testthat::test_that("validate_lc_column_id function works correctly", {
  
  # Test data
  df <- data.frame(
    lc_column_id = c("LC01", "LC02", "LC03", "LC 04", "LC05", NA, ""),
    other_column = 1:7
  )
  
  # Verify outcomes
  testthat::expect_equal(MotrpacBicQC::validate_lc_column_id(df, "lc_column_id", verbose = FALSE), 3)
  
})

testthat::test_that("validate_lc_column_id function handles missing column", {
  
  # Test data
  df <- data.frame(other_column = 1:3)
  
  # Run function under test and verify error
  testthat::expect_error(MotrpacBicQC::validate_lc_column_id(df, "lc_column_id", verbose = TRUE),
                         "Column lc_column_id does not exist in the data frame.")
  
})

testthat::test_that("validate_lc_column_id function handles excess unique values", {
  
  # Test data
  df <- data.frame(
    lc_column_id = c("LC01", "LC02", "LC03", "LC04", "LC05", "LC06", "LC07"),
    other_column = 1:7
  )
  
  # Verify outcomes
  testthat::expect_equal(MotrpacBicQC::validate_lc_column_id(df, "lc_column_id", verbose = FALSE), 0)
  
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test case 1: The column exists and there are no NA or empty values
testthat::test_that("validate_na_empty function handles no NA or empty values", {
  df <- data.frame(
    existing_column = c("value1", "value2", "value3")
  )
  testthat::expect_equal(MotrpacBicQC::validate_na_empty(df, "existing_column", verbose = FALSE), 0)
})

# Test case 2: The column exists and there are NA values
testthat::test_that("validate_na_empty function handles NA values", {
  df <- data.frame(
    existing_column = c("value1", NA, "value3")
  )
  testthat::expect_equal(MotrpacBicQC::validate_na_empty(df, "existing_column", verbose = FALSE), 1)
})

# Test case 3: The column exists and there are empty values
testthat::test_that("validate_na_empty function handles empty values", {
  df <- data.frame(
    existing_column = c("value1", "", "value3")
  )
  testthat::expect_equal(MotrpacBicQC::validate_na_empty(df, "existing_column", verbose = FALSE), 1)
})

# Test case 4: The column does not exist
testthat::test_that("validate_na_empty function handles non-existing column", {
  df <- data.frame(
    existing_column = c("value1", "value2", "value3")
  )
  testthat::expect_error(MotrpacBicQC::validate_na_empty(df, "non_existing_column", verbose = FALSE), 
                         "This column non_existing_column does not exist")
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test case 1: The column exists and all dates are in the correct format and range
testthat::test_that("validate_yyyymmdd_dates function handles correct dates", {
  df <- data.frame(
    existing_column = c("2023-01-01", "2023-12-31", "2024-02-29")
  )
  testthat::expect_equal(MotrpacBicQC::validate_yyyymmdd_dates(df, "existing_column", verbose = FALSE), 0)
})

# Test case 2: The column exists and contains dates with incorrect format
testthat::test_that("validate_yyyymmdd_dates function handles incorrect format", {
  df <- data.frame(
    existing_column = c("2023-01-01", "20231231", "2024-02-29")
  )
  testthat::expect_equal(MotrpacBicQC::validate_yyyymmdd_dates(df, "existing_column", verbose = FALSE), 1)
})

# Test case 3: The column exists and contains dates with incorrect components
testthat::test_that("validate_yyyymmdd_dates function handles incorrect components", {
  df <- data.frame(
    existing_column = c("2023-01-01", "2023-13-31", "2024-02-30")
  )
  testthat::expect_equal(MotrpacBicQC::validate_yyyymmdd_dates(df, "existing_column", verbose = FALSE), 1)
})

# Test case 4: The column does not exist
testthat::test_that("validate_yyyymmdd_dates function handles non-existing column", {
  df <- data.frame(
    existing_column = c("2023-01-01", "2023-12-31", "2024-02-29")
  )
  testthat::expect_error(MotrpacBicQC::validate_yyyymmdd_dates(df, "non_existing_column", verbose = FALSE), 
                         "Column non_existing_column not found in the data frame.")
})

# Test case 5: The column exists and contains dates with / instead of -
testthat::test_that("validate_yyyymmdd_dates function handles dates with / instead of -", {
  df <- data.frame(
    existing_column = c("2023/01/01", "2023-12-31", "2024-02-29")
  )
  testthat::expect_equal(MotrpacBicQC::validate_yyyymmdd_dates(df, "existing_column", verbose = FALSE), 1)
})



