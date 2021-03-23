context("Test Validations")

test_that("All validations works", {
  x <- "PASS1A-06/T31/IONPNEG/BATCH1_20190909/PROCESSED_20200205"
  y <- "MORESTUFF/PASS1A-06/T31/IONPNEG/BATCH1_20190909/PROCESSED_20200205/EVENMOREFOLDERS/"
  z <- "PASS1B-06/T100/IONPNEG/BATCH1_20190909/PROCESSED_20200205"
  expect_equal(validate_processFolder(x), "PROCESSED_20200205")
  expect_equal(validate_processFolder(y), "PROCESSED_20200205")
  expect_equal(validate_assay(x), "IONPNEG")
  expect_equal(validate_phase(x), "PASS1A-06")
  expect_equal(validate_phase(y), "PASS1A-06")
  expect_equal(validate_tissue(x), "T31")
  expect_error(validate_tissue(z))
})


