context("Test Validations")

test_that("All validations works", {
  x <- "PASS1A-06/T31/IONPNEG/BATCH1_20190909/PROCESSED_20200205"
  y <- "MORESTUFF/PASS1A-06/T31/IONPNEG/BATCH1_20190909/PROCESSED_20200205/EVENMOREFOLDERS/"
  z <- "PASS1B-06/T100/IONPNEG/BATCH1_20190909/PROCESSED_20200205"
  b <- "PASS1A-06/T31/IONPNEG/BATCH1_2019090990/PROCESSED_20200205"
  h <- "HUMAN/T10/IONPNEG/BATCH1_20190909/RESULTS_20220101"
  j <- "HUMAN/T10/IONPNEG/BATCH1_20190909/BICRESULTS_20220101"
  k <- "HUMAN/T10/IONPNEG/BATCH1_20190909/TOPRESULTS_20220101"
  l <- "HUMAN/PRECOVID/T10/PROT_PR/BATCH1_20220919/PROCESSED_20200205"
  m <- "PASS1A-06|PASS1C-06"
  expect_equal(validate_processFolder(x), "PROCESSED_20200205")
  expect_equal(validate_processFolder(y), "PROCESSED_20200205")
  expect_equal(validate_assay(x), "IONPNEG")
  expect_equal(validate_phase(x), "PASS1A-06")
  expect_equal(validate_phase(y), "PASS1A-06")
  expect_equal(validate_tissue(x), "T31")
  expect_equal(validate_batch(h), "HUMAN/T10/IONPNEG/BATCH1_20190909/")
  expect_equal(validate_batch(j), "HUMAN/T10/IONPNEG/BATCH1_20190909/")
  expect_equal(validate_processFolder(h), "RESULTS_20220101")
  expect_equal(validate_processFolder(j), "BICRESULTS_20220101")
  expect_error(expect_match(validate_processFolder(k), "TOPRESULTS_20220101"))
  expect_error(validate_tissue(z))
  expect_equal(validate_batch(x), "PASS1A-06/T31/IONPNEG/BATCH1_20190909/")
  expect_error(validate_batch(b))
  expect_equal(generate_phase_details(phase_metadata = m), "pass1ac-06")
  expect_error(set_phase(m))
  expect_equal(validate_two_phases(m, verbose = TRUE), "Two phases reported and they are ok")
  expect_error(validate_two_phases(phase_details = "PASS1A-06"))
})


