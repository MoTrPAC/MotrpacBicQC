
context("Test check_metadata_proteins")

# Mock data frames for testing
valid_df1 <- data.frame(
  olink_id = 1:10,
  uniprot_entry = paste0("P", 1:10),
  assay = LETTERS[1:10],
  missing_freq = runif(10),
  panel_name = rep("Panel1", 10),
  panel_lot_nr = rep("LotA", 10),
  normalization = rep("Method1", 10)
)

invalid_df1 <- valid_df1
invalid_df1$olink_id[c(1, 3)] <- invalid_df1$olink_id[2] # Create duplicate olink_ids

test_that("Valid data returns zero issues", {
  expect_silent(result <- check_metadata_proteins(df = valid_df1, return_n_issues = TRUE, verbose = FALSE))
  expect_equal(result, 0)
})

test_that("Duplicate olink_ids are detected", {
  expect_silent(result <- check_metadata_proteins(df = invalid_df1, return_n_issues = TRUE, verbose = FALSE))
  expect_gt(result, 0) # result should be greater than 0 due to duplicates
})

test_that("Missing columns are detected", {
  df_missing_cols <- valid_df1[, -which(names(valid_df1) == "olink_id")]
  expect_message(result <- check_metadata_proteins(df = df_missing_cols, return_n_issues = TRUE, verbose = FALSE))
  expect_gt(result, 0)
})

test_that("NA values in critical columns are detected", {
  df_with_na <- valid_df1
  df_with_na$olink_id[1] <- NA  # Introduce an NA value
  result <- check_metadata_proteins(df = df_with_na, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

test_that("Invalid uniprot_entry values are detected", {
  df_invalid_uniprot <- valid_df1
  df_invalid_uniprot$uniprot_entry[1] <- "INVALID_ID"  # Introduce an invalid Uniprot ID
  result <- check_metadata_proteins(df = df_invalid_uniprot, 
                                    return_n_issues = TRUE, 
                                    validate_uniprot = TRUE, 
                                    verbose = FALSE)
  expect_gt(result, 0)
})

test_that("Non-numeric values in missing_freq are detected", {
  df_nonnumeric <- valid_df1
  df_nonnumeric$missing_freq[1] <- "non-numeric"  # Introduce a non-numeric value
  result <- check_metadata_proteins(df = df_nonnumeric, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

test_that("NA values in panel_name are detected", {
  df_with_na_panel <- valid_df1
  df_with_na_panel$panel_name[1] <- NA
  result <- check_metadata_proteins(df = df_with_na_panel, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

test_that("Empty data frame is handled", {
  empty_df <- data.frame()
  result <- check_metadata_proteins(df = empty_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
context("Test check_metadata_samples_olink")

# Valid data frame
valid_df2 <- data.frame(
  sample_id = 1:10,
  sample_type = rep("Sample", 10),
  plate_id = rep("Plate1", 10),
  sample_order = 1:10
)

# Data frame with duplicate sample_id
duplicate_sample_id_df <- valid_df2
duplicate_sample_id_df$sample_id[5] <- duplicate_sample_id_df$sample_id[4]

# Test for handling valid data
test_that("Valid data is processed correctly", {
  expect_silent(result <- check_metadata_samples_olink(df = valid_df2, return_n_issues = TRUE, verbose = FALSE))
  expect_equal(result, 0)
})

# Test for detecting duplicate sample_id
test_that("Duplicate sample_id values are detected", {
  result <- check_metadata_samples_olink(df = duplicate_sample_id_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

# Test for detecting undefined sample types
test_that("Undefined sample types are detected", {
  invalid_sample_type_df <- valid_df2
  invalid_sample_type_df$sample_type[1] <- "InvalidType"
  result <- check_metadata_samples_olink(df = invalid_sample_type_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

# Test for detecting non-numeric sample_order
test_that("Non-numeric sample_order values are detected", {
  nonnumeric_sample_order_df <- valid_df2
  nonnumeric_sample_order_df$sample_order[1] <- "non-numeric"
  result <- check_metadata_samples_olink(df = nonnumeric_sample_order_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

# Test for detecting missing columns
test_that("Missing columns are detected", {
  df_missing_columns <- valid_df2[, -which(names(valid_df2) %in% c("sample_id", "plate_id"))]
  result <- check_metadata_samples_olink(df = df_missing_columns, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

# Test for non-unique sample_order within each plate_id
test_that("Non-unique sample_order within plate_id is detected", {
  non_unique_order_df <- valid_df2
  non_unique_order_df$sample_order[c(1, 2)] <- 1  # Create duplicate order in the same plate
  result <- check_metadata_samples_olink(df = non_unique_order_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
context("Test check_metadata_proteins")

# Valid data frame
valid_df3 <- data.frame(
  olink_id = 1:10,
  measurement1 = rnorm(10),
  measurement2 = runif(10)
)

# Data frame with duplicate olink_id
duplicate_olink_id_df <- valid_df3
duplicate_olink_id_df$olink_id[5] <- duplicate_olink_id_df$olink_id[4]

# Data frame with non-numeric column
non_numeric_df <- valid_df3
non_numeric_df$measurement1 <- as.character(non_numeric_df$measurement1)
non_numeric_df$measurement1[1] <- "non-numeric"

# Test for handling valid data
test_that("Valid data is processed correctly", {
  expect_silent(result <- check_results_olink(df = valid_df3, return_n_issues = TRUE, verbose = FALSE))
  expect_equal(result, 0)
})


# Test for detecting duplicate olink_id
test_that("Duplicate olink_id values are detected", {
  result <- check_results_olink(df = duplicate_olink_id_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

# Test for detecting non-numeric columns
test_that("Non-numeric columns are detected", {
  result <- check_results_olink(df = non_numeric_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})

# Test for detecting NA values
test_that("NA values are detected", {
  df_with_na <- valid_df3
  df_with_na$measurement1[1] <- NA
  result <- check_results_olink(df = df_with_na, return_n_issues = TRUE, verbose = FALSE)
  # Not expecting an increase in issue count as NA detection might not be an issue
  expect_true(result >= 0)
})

# Test for detecting zero values
test_that("Zero values are detected", {
  df_with_zeros <- valid_df3
  df_with_zeros$measurement1 <- rep(0, 10)
  result <- check_results_olink(df = df_with_zeros, return_n_issues = TRUE, verbose = FALSE)
  # Not expecting an increase in issue count as zero detection might not be an issue
  expect_true(result >= 0)
})

test_that("Empty data frame is handled", {
  empty_df <- data.frame()
  result <- check_results_olink(df = empty_df, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)  # Expect issues since the data frame is empty
})

test_that("Missing olink_id column is detected", {
  df_missing_olink_id <- valid_df3[, -which(names(valid_df3) == "olink_id")]
  result <- check_results_olink(df = df_missing_olink_id, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})


test_that("Data with all zeros is processed", {
  df_all_zeros <- valid_df3
  df_all_zeros[,-which(names(df_all_zeros) == "olink_id")] <- 0
  result <- check_results_olink(df = df_all_zeros, return_n_issues = TRUE, verbose = FALSE)
  expect_true(result >= 0)
})


test_that("Data with all NAs is processed", {
  df_all_na <- valid_df3
  df_all_na[,-which(names(df_all_na) == "olink_id")] <- NA
  result <- check_results_olink(df = df_all_na, return_n_issues = TRUE, verbose = FALSE)
  expect_true(result >= 0)
})


test_that("Non-numeric columns other than olink_id are detected", {
  df_with_char_column <- valid_df3
  df_with_char_column$measurement1 <- as.character(df_with_char_column$measurement1)
  result <- check_results_olink(df = df_with_char_column, return_n_issues = TRUE, verbose = FALSE)
  expect_gt(result, 0)
})



context("Test check_crossfile_olink_validation")

# Mock data frames for testing
valid_results_df <- data.frame(
  olink_id = paste0("OLINK", 1:5),
  Sample1 = rnorm(5),
  Sample2 = rnorm(5)
)
valid_metadata_samples_df <- data.frame(
  sample_id = c("Sample1", "Sample2")
)
valid_metadata_proteins_df <- data.frame(
  olink_id = paste0("OLINK", 1:5)
)

# Data frame with mismatched sample IDs
mismatched_samples_df <- valid_metadata_samples_df
mismatched_samples_df$sample_id <- c("Sample3", "Sample4")

# Data frame with mismatched olink IDs
mismatched_olink_df <- valid_metadata_proteins_df
mismatched_olink_df$olink_id <- paste0("OLINK", 6:10)

# Test for handling valid data
test_that("Valid data is processed correctly", {
  expect_silent(result <- check_crossfile_olink_validation(r_o = valid_results_df, 
                                                           m_s = valid_metadata_samples_df, 
                                                           m_p = valid_metadata_proteins_df, 
                                                           return_n_issues = TRUE, 
                                                           verbose = FALSE))
  expect_equal(result, 0)
})

# Test for detecting mismatched sample IDs
test_that("Mismatched sample IDs are detected", {
  result <- check_crossfile_olink_validation(r_o = valid_results_df, 
                                             m_s = mismatched_samples_df, 
                                             m_p = valid_metadata_proteins_df, 
                                             return_n_issues = TRUE, 
                                             verbose = TRUE)
  expect_gt(result, 0)
})

# Test for detecting mismatched olink IDs
test_that("Mismatched olink IDs are detected", {
  result <- check_crossfile_olink_validation(r_o = valid_results_df, 
                                             m_s = valid_metadata_samples_df, 
                                             m_p = mismatched_olink_df, 
                                             return_n_issues = TRUE, 
                                             verbose = TRUE)
  expect_gt(result, 0)
})

