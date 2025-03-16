library(testthat)
library(MotrpacBicQC)

context("QC and Data Processing Functions")

# --- Dummy Data Setup ---

dummy_analyte <- data.frame(
  analyte_name   = c("A1", "A2", "A3", "A1"),  # duplicate "A1"
  uniprot_entry  = c("P12345", "Q67890", "A11111", "P12345"),
  assay_name     = rep("assay1", 4),
  stringsAsFactors = FALSE
)

dummy_analyte_unique <- data.frame(
  analyte_name   = c("A1", "A2", "A3"),
  uniprot_entry  = c("P12345", "Q67890", "A11111"),
  assay_name     = rep("assay1", 3),
  stringsAsFactors = FALSE
)

dummy_samples <- data.frame(
  sample_id      = c("S1", "S2", "S3", "S3"),  # duplicate "S3"
  sample_type    = c("Sample", "QC-Pooled", "QC-Pooled", "Sample"),
  sample_order   = c(1, 2, 3, 4),
  raw_file       = c("file1.txt", "file2.txt", "file3.txt", "file4.txt"),
  extraction_date= c("20200101", "20200102", "20200103", "20200104"),
  acquisition_date = c("20200105", "20200106", "20200107", "20200108"),
  stringsAsFactors = FALSE
)

# --- Tests for check_metadata_analyte ---

test_that("check_metadata_analyte detects duplicate analyte_name", {
  issues <- check_metadata_analyte(dummy_analyte, return_n_issues = TRUE, verbose = FALSE)
  expect_true(issues >= 1)
})

test_that("check_metadata_analyte passes with unique analyte_name", {
  issues <- check_metadata_analyte(dummy_analyte_unique, return_n_issues = TRUE, verbose = FALSE)
  expect_equal(issues, 0)
})

# --- Tests for check_metadata_samples_lab ---

test_that("check_metadata_samples_lab detects non-unique sample_id", {
  issues <- check_metadata_samples_lab(dummy_samples, return_n_issues = TRUE, verbose = FALSE)
  expect_true(issues >= 1)
})

test_that("check_metadata_analyte validates analyte metadata correctly", {
  # Create a valid metadata_analyte dataframe
  valid_df <- data.frame(
    analyte_name = c("Analyte1", "Analyte2", "Analyte3"),
    uniprot_entry = c("P12345", "P23456", "P34567"),
    assay_name = c("Assay1", "Assay2", "Assay3"),
    stringsAsFactors = FALSE
  )

  # Test with valid data
  expect_equal(check_metadata_analyte(valid_df, return_n_issues = TRUE, verbose = FALSE), 0)

  # Test with missing analyte_name column
  invalid_df <- valid_df[, c("uniprot_entry", "assay_name")]
  expect_gt(check_metadata_analyte(invalid_df, return_n_issues = TRUE, verbose = FALSE), 0)

  # Test with duplicate analyte_name
  duplicate_df <- rbind(valid_df, valid_df[1,])
  expect_gt(check_metadata_analyte(duplicate_df, return_n_issues = TRUE, verbose = FALSE), 0)

  # Test with NA values
  na_df <- valid_df
  na_df$analyte_name[1] <- NA
  expect_gt(check_metadata_analyte(na_df, return_n_issues = TRUE, verbose = FALSE), 0)
})

test_that("check_metadata_samples_lab validates sample metadata correctly", {
  # Create a valid metadata_sample dataframe
  valid_df <- data.frame(
    sample_id = c("Sample1", "Sample2", "Sample3", "QC1"),
    sample_type = c("Sample", "Sample", "Sample", "QC-Pooled"),
    sample_order = c(1, 2, 3, 4),
    raw_file = c("file1.raw", "file2.raw", "file3.raw", "file4.raw"),
    extraction_date = c("2023-01-01", "2023-01-01", "2023-01-01", "2023-01-01"),
    acquisition_date = c("2023-01-02", "2023-01-02", "2023-01-02", "2023-01-02"),
    stringsAsFactors = FALSE
  )

  # Test with valid data
  expect_equal(check_metadata_samples_lab(valid_df, return_n_issues = TRUE, verbose = FALSE), 0)

  # Test with missing sample_id column
  invalid_df <- valid_df[, !names(valid_df) %in% "sample_id"]
  expect_gt(check_metadata_samples_lab(invalid_df, return_n_issues = TRUE, verbose = FALSE), 0)

  # Test with invalid sample_type
  invalid_type_df <- valid_df
  invalid_type_df$sample_type[1] <- "InvalidType"
  expect_gt(check_metadata_samples_lab(invalid_type_df, return_n_issues = TRUE, verbose = FALSE), 0)

  # Test with non-numeric sample_order
  non_numeric_df <- valid_df
  non_numeric_df$sample_order[1] <- "a"
  expect_gt(check_metadata_samples_lab(non_numeric_df, return_n_issues = TRUE, verbose = FALSE), 0)

  # Test with NA in raw_file
  na_raw_file_df <- valid_df
  na_raw_file_df$raw_file[1] <- NA
  expect_gt(check_metadata_samples_lab(na_raw_file_df, return_n_issues = TRUE, verbose = FALSE), 0)
})

test_that("load_lab_batch loads data correctly", {
  # Mock the required functions using local_mocked_bindings instead of with_mock
  local_mocked_bindings(
    validate_phase = function(...) "PASS1B",
    validate_processFolder = function(...) "PROCESSED_20230101",
    validate_assay = function(...) "CK",
    validate_tissue = function(...) "T01",
    validate_lab = function(...) 0,
    open_file = function(input_results_folder, filepattern, verbose) {
      if(grepl("metadata_analyte", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "metadata_analyte_test.txt",
          df = data.frame(
            analyte_name = c("Analyte1", "Analyte2"),
            uniprot_entry = c("P12345", "P23456"),
            assay_name = c("Assay1", "Assay2")
          )
        ))
      } else if(grepl("metadata_sample", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "metadata_sample_test.txt",
          df = data.frame(
            sample_id = c("Sample1", "Sample2"),
            sample_type = c("Sample", "QC-Pooled"),
            sample_order = c(1, 2),
            raw_file = c("file1.raw", "file2.raw"),
            extraction_date = as.Date(c("2023-01-01", "2023-01-01")),
            acquisition_date = as.Date(c("2023-01-02", "2023-01-02"))
          )
        ))
      } else if(grepl("results", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "results_test.txt",
          df = data.frame(
            analyte_name = c("Analyte1", "Analyte2"),
            Sample1 = c(10.5, 11.2),
            Sample2 = c(12.3, 13.1)
          )
        ))
      }
    }
  )

  result <- load_lab_batch(input_results_folder = "mock/path", verbose = FALSE)

  expect_type(result, "list")
  expect_equal(names(result), c("m_a", "m_s", "r_o"))
  expect_equal(nrow(result$m_a), 2)
  expect_equal(nrow(result$m_s), 2)
  expect_equal(nrow(result$r_o), 2)
})

test_that("validate_lab correctly identifies issues", {
  skip_if_not_installed("mockery")

  local_mocked_bindings(
    validate_cas = function(...) NULL,
    validate_processFolder = function(...) "PROCESSED_20230101",
    validate_assay = function(...) "CK",
    validate_phase = function(...) "PASS1B",
    validate_tissue = function(...) "T01",
    validate_batch = function(...) "BATCH1_20230101",
    check_metadata_phase_file = function(...) TRUE,
    set_phase = function(...) "PASS1B",
    generate_phase_details = function(...) "pass1b",
    open_file = function(input_results_folder, filepattern, verbose) {
      if(grepl("metadata_analyte", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "metadata_analyte_test.txt",
          df = data.frame(
            analyte_name = c("Analyte1", "Analyte2"),
            uniprot_entry = c("P12345", "P23456"),
            assay_name = c("Assay1", "Assay2")
          )
        ))
      } else if(grepl("metadata_sample", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "metadata_sample_test.txt",
          df = data.frame(
            sample_id = c("Sample1", "Sample2"),
            sample_type = c("Sample", "QC-Pooled"),
            sample_order = c(1, 2),
            raw_file = c("file1.raw", "file2.raw"),
            extraction_date = as.Date(c("2023-01-01", "2023-01-01")),
            acquisition_date = as.Date(c("2023-01-02", "2023-01-02"))
          )
        ))
      } else if(grepl("results", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "results_test.txt",
          df = data.frame(
            analyte_name = c("Analyte1", "Analyte2"),
            Sample1 = c(10.5, 11.2),
            Sample2 = c(12.3, 13.1)
          )
        ))
      }
    },
    check_metadata_analyte = function(...) 0,
    check_metadata_samples_lab = function(...) 0,
    check_results_assays = function(...) 0,
    check_crossfile_validation = function(...) 0,
    check_failedsamples = function(...) NULL
  )

  mockery::stub(validate_lab, "normalizePath", function(...) "/mock/path/HUMAN/T01/LAB_CK/BATCH1_20230101/PROCESSED_20230101")
  mockery::stub(validate_lab, "regmatches", function(...) "HUMAN/T01/LAB_CK/BATCH1_20230101/PROCESSED_20230101")
  mockery::stub(validate_lab, "list.files", function(...) character(0))
  mockery::stub(validate_lab, "purrr::is_empty", function(...) FALSE)

  # Run the test
  result <- validate_lab(
    input_results_folder = "mock/path",
    cas = "duke",
    return_n_issues = TRUE,
    verbose = FALSE
  )

  # Check the result for success case
  expect_type(result, "double")
  expect_lte(result, 10)  # Should report few or no issues
})

test_that("validate_lab identifies issues with missing files", {
  
  skip_if_not_installed("mockery")

  # Mock the required functions
  local_mocked_bindings(
    validate_cas = function(...) NULL,
    validate_processFolder = function(...) "PROCESSED_20230101",
    validate_assay = function(...) "CK",
    validate_phase = function(...) "PASS1B",
    validate_tissue = function(...) "T01",
    validate_batch = function(...) "BATCH1_20230101",
    check_metadata_phase_file = function(...) TRUE,
    set_phase = function(...) "PASS1B",
    generate_phase_details = function(...) "pass1b",
    open_file = function(input_results_folder, filepattern, verbose) {
      if(grepl("metadata_analyte", filepattern)) {
        return(list(flag = FALSE))  # Missing file
      } else if(grepl("metadata_sample", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "metadata_sample_test.txt",
          df = data.frame(
            sample_id = c("Sample1", "Sample2"),
            sample_type = c("Sample", "QC-Pooled"),
            sample_order = c(1, 2),
            raw_file = c("file1.raw", "file2.raw"),
            extraction_date = as.Date(c("2023-01-01", "2023-01-01")),
            acquisition_date = as.Date(c("2023-01-02", "2023-01-02"))
          )
        ))
      } else if(grepl("results", filepattern)) {
        return(list(
          flag = TRUE,
          filename = "results_test.txt",
          df = data.frame(
            analyte_name = c("Analyte1", "Analyte2"),
            Sample1 = c(10.5, 11.2),
            Sample2 = c(12.3, 13.1)
          )
        ))
      }
    },
    check_metadata_samples_lab = function(...) 0,
    check_results_assays = function(...) 0,
    check_failedsamples = function(...) NULL
  )

  # Use mockery for base R and external package functions
  mockery::stub(validate_lab, "normalizePath", function(...) "/mock/path/HUMAN/T01/LAB_CK/BATCH1_20230101/PROCESSED_20230101")
  mockery::stub(validate_lab, "regmatches", function(...) "HUMAN/T01/LAB_CK/BATCH1_20230101/PROCESSED_20230101")
  mockery::stub(validate_lab, "list.files", function(...) character(0))
  mockery::stub(validate_lab, "purrr::is_empty", function(...) FALSE)

  # Run the test
  result <- validate_lab(
    input_results_folder = "mock/path",
    cas = "duke",
    return_n_issues = TRUE,
    verbose = FALSE
  )

  # Check the result for failure case
  expect_type(result, "double")
  expect_gt(result, 10)  # Should report significant issues
})

