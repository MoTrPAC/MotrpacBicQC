context("Test misc functions")

test_that("remove empty columns works", {
  L3 <- LETTERS[1:3]
  fac <- sample(L3, 10, replace = TRUE)
  df <- data.frame(x = 1, y = 1:10, fac = fac, stringsAsFactors = FALSE)

  # Add empty column
  df$tremove1 <- ""
  df$tremove2 <- NA

  expect_equal(dim(remove_empty_columns(df, verbose = FALSE))[2], 3)
})

test_that("remove empty rows works", {
  L3 <- LETTERS[1:3]
  fac <- sample(L3, 10, replace = TRUE)
  df <- data.frame(x = 1, y = 1:10, fac = fac, stringsAsFactors = FALSE)

  # Add empty rows
  df[nrow(df)+1,] <- NA
  df[nrow(df)+1,] <- ""

  expect_equal(dim(remove_empty_rows(df, verbose = FALSE))[1], 10)
})

test_that("Check all columns in data table bic_animal_tissue_code", {
  expected_colnames <- c("bic_tissue_code", "bic_tissue_name", "motrpac_tissue_code", "tissue_name_release", "abbreviation", "tissue_hex_colour")
  expect_equal( setequal(expected_colnames, colnames(bic_animal_tissue_code)), TRUE)
  expect_gt(dim(bic_animal_tissue_code)[1], 40)
})

test_that("Whitespace is trimmed in character columns", {
  df <- data.frame(
    a = c("  hello", "world  ", "  foo  "),
    b = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  cleaned <- clean_character_columns(df)
  expect_equal(cleaned$a, c("hello", "world", "foo"))
  expect_equal(cleaned$b, c(1, 2, 3))
})

test_that("Function leaves columns without extra whitespace unchanged", {
  df <- data.frame(
    a = c("apple", "banana", "cherry"),
    b = c("dog", "cat", "bird"),
    stringsAsFactors = FALSE
  )
  cleaned <- clean_character_columns(df)
  expect_equal(cleaned$a, c("apple", "banana", "cherry"))
  expect_equal(cleaned$b, c("dog", "cat", "bird"))
})

test_that("Works with empty data frame", {
  df <- data.frame(
    a = character(0),
    b = numeric(0),
    stringsAsFactors = FALSE
  )
  cleaned <- clean_character_columns(df)
  expect_equal(nrow(cleaned), 0)
})

