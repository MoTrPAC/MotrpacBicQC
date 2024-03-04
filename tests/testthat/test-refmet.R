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


test_that("Successful API call returns correctly structured list with no additional elements", {
  # Known good refmet_name that will return a successful response
  test_refmet_name <- "11-Deoxycortisol"
  
  # Expected elements in the response
  expected_elements <- c("refmet_name", "formula", "exactmass", "super_class", "main_class", "sub_class")
  
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

