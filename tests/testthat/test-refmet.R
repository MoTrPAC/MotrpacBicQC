context("Testing validate_refmetname function")

test_that("validate_refmetname handles known refmet_name correctly", {
  # Example test data
  test_data <- data.frame(
    refmet_name = c("1-Methyl nicotinamide",  "11-Deoxycortisol", "This Should Fail", "I hope this as well---", "-"),
    stringsAsFactors = FALSE
  )
  
  # Call the function with verbose = FALSE to suppress messages during test
  actual_missed_ids <- validate_refmetname(test_data, verbose = FALSE)
  
  # Check if the actual missed IDs match the expected outcome
  expect_equal(actual_missed_ids, 3)
})

test_that("validate_refmetname to very similar refmet_names", {
  # Example test data
  # https://www.metabolomicsworkbench.org/rest/refmet/match/2-Amino-6-hydroxyhexanoic acid/name/
  # https://www.metabolomicsworkbench.org/rest/refmet/match/2-amino-6-hydroxyhexanoic acid/name/
  test_data <- data.frame(
    refmet_name = c("2-Amino-6-hydroxyhexanoic acid",  "2-amino-6-hydroxyhexanoic acid", "2-Aamino-6-hydroxyhexanoic acid", "-"),
    stringsAsFactors = FALSE
  )
  
  # Call the function with verbose = FALSE to suppress messages during test
  actual_missed_ids <- validate_refmetname(test_data, verbose = TRUE)
  
  # Check if the actual missed IDs match the expected outcome
  expect_equal(actual_missed_ids, 3)
})


test_that("Successful API call returns correctly structured list with no additional elements", {
  # Known good refmet_name that will return a successful response
  test_refmet_name <- "11-Deoxycortisol"
  
  # Expected elements in the response
  expected_elements <- c("refmet_name", "formula", "exactmass", "super_class", "main_class", "sub_class", "refmet_id")
  
  search_api <- paste0("https://www.metabolomicsworkbench.org/rest/refmet/match/",URLencode(test_refmet_name),"/name/")
  response <- jsonlite::fromJSON(search_api)
  
  # Check that the response has all the expected elements
  expect_true(all(expected_elements %in% names(response)), "Response is missing expected elements.")
  
  # Check for no additional elements
  expect_equal(length(names(response)), length(expected_elements), 
               info = "Response contains additional unexpected elements.")
  
  unexpected_elements <- setdiff(names(response), expected_elements)
  expect_equal(unexpected_elements, character(0), 
               info = paste("Unexpected elements in response:", paste(unexpected_elements, collapse = ", ")))
})

test_that("validate_refmetname accepts known slash-containing names", {
  skip_if_offline()

  test_data <- data.frame(
    refmet_name = c("Leucine/Isoleucine", "Citric acid/Isocitric acid"),
    stringsAsFactors = FALSE
  )

  actual_missed_ids <- validate_refmetname(test_data, verbose = FALSE)
  expect_equal(actual_missed_ids, 0)
})

test_that("validate_refmetname skips [iSTD] entries", {
  test_data <- data.frame(
    refmet_name = c("alanine-13C3 [iSTD]", "anthranilic acid-15N [iSTD]", "24-epi-Brassinolide [iSTD]"),
    stringsAsFactors = FALSE
  )
  
  actual_missed_ids <- validate_refmetname(test_data, verbose = FALSE)
  expect_equal(actual_missed_ids, 0)
})

test_that("validate_refmetname rejects non-data.frame inputs instead of silently passing", {
  # Passing a plain list used to silently return 0 because nrow(list) is NULL.
  # A list carrying a `refmet_name` element should be coerced and validated.
  skip_if_offline()
  expect_equal(
    validate_refmetname(list(refmet_name = c("Leucine", "Citric acid")), verbose = FALSE),
    0
  )
  # Inputs without a usable `refmet_name` column must error loudly.
  expect_error(
    validate_refmetname(c("Leucine", "Citric acid"), verbose = FALSE),
    "data.frame"
  )
  expect_error(
    validate_refmetname(list(other = c("Leucine", "Citric acid")), verbose = FALSE),
    "data.frame"
  )
  expect_error(
    validate_refmetname(data.frame(other = "x", stringsAsFactors = FALSE), verbose = FALSE),
    "refmet_name"
  )
})

test_that("validate_refmetname returns 0 for an empty data.frame", {
  empty_df <- data.frame(refmet_name = character(0), stringsAsFactors = FALSE)
  expect_equal(validate_refmetname(empty_df, verbose = FALSE), 0)
})

test_that("validate_refmetname does not crash on API/network errors", {
  skip_if_not_installed("mockery")

  # Simulate jsonlite::fromJSON failing for every call (network/API error).
  mockery::stub(
    validate_refmetname,
    "jsonlite::fromJSON",
    function(...) stop("simulated network failure")
  )

  test_data <- data.frame(
    refmet_name = c("Leucine", "Citric acid"),
    stringsAsFactors = FALSE
  )

  result <- validate_refmetname(test_data, verbose = FALSE)
  # Both rows could not be verified, so they are counted as missed IDs
  # rather than aborting the whole QC run.
  expect_equal(result, 2)
})

test_that("validate_refmetname handles NA/empty refmet_name without crashing", {
  # Mock the API so the test doesn't depend on the network. Only non-NA/
  # non-empty rows should reach the API; missing rows must be counted
  # locally as Error RN0.
  skip_if_not_installed("mockery")
  mockery::stub(
    validate_refmetname,
    "jsonlite::fromJSON",
    function(url, ...) list(refmet_name = "Leucine")
  )

  test_data <- data.frame(
    refmet_name = c("Leucine", NA_character_, "", "   "),
    stringsAsFactors = FALSE
  )

  # 3 missing/empty rows => counted as missed IDs (Error RN0), but no crash.
  result <- validate_refmetname(test_data, verbose = FALSE)
  expect_equal(result, 3)
})

