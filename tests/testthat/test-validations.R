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
  # PROT_AC legacy exception: BATCH_YYYYMMDD (no batch number) is allowed
  prot_ac_legacy <- "broad/PASS1A-06/T58/PROT_AC/BATCH_20190828/PROCESSED_20200101"
  expect_equal(validate_batch(prot_ac_legacy), "broad/PASS1A-06/T58/PROT_AC/BATCH_20190828/")
  # PROT_AC with a batch number should still work
  prot_ac_numbered <- "broad/PASS1A-06/T58/PROT_AC/BATCH1_20190828/PROCESSED_20200101"
  expect_equal(validate_batch(prot_ac_numbered), "broad/PASS1A-06/T58/PROT_AC/BATCH1_20190828/")
  # Non-PROT_AC without batch number must fail
  non_prot_no_num <- "broad/PASS1A-06/T58/IONPNEG/BATCH_20190828/PROCESSED_20200101"
  expect_error(validate_batch(non_prot_no_num))
  expect_equal(generate_phase_details(phase_metadata = m), "pass1ac-06")
  expect_error(set_phase(m))
  expect_equal(validate_two_phases(m, verbose = TRUE), "Two phases reported and they are ok")
  expect_error(validate_two_phases(phase_details = "PASS1A-06"))
})

test_that("validate_phase rejects HUMAN-MAIN-TR## in folder path", {
  # HUMAN-MAIN-TR01 as a folder should fail (must be just HUMAN)
  expect_error(validate_phase("HUMAN-MAIN-TR01/T11/OXYLIPNEG/BATCH1_20240524/PROCESSED_20240524"))
  # HUMAN-MAIN (without tranche) should also fail

  expect_error(validate_phase("HUMAN-MAIN/T11/OXYLIPNEG/BATCH1_20240524/PROCESSED_20240524"))
  # HUMAN_MAIN_TR01 (underscores) should fail
  expect_error(validate_phase("HUMAN_MAIN_TR01/T11/OXYLIPNEG/BATCH1_20240524/PROCESSED_20240524"))
})

test_that("validate_phase accepts valid folder phases", {
  expect_equal(validate_phase("HUMAN/T11/OXYLIPNEG/BATCH1_20240524/PROCESSED_20240524"), "HUMAN")
  expect_equal(validate_phase("HUMAN-PRECOVID/T10/IONPNEG/BATCH1_20190909/PROCESSED_20200205"), "HUMAN-PRECOVID")
  expect_equal(validate_phase("PASS1A-06/T31/IONPNEG/BATCH1_20190909/PROCESSED_20200205"), "PASS1A-06")
  # Bare phase strings (used by validate_two_phases) should also work
  expect_equal(validate_phase("PASS1C-06"), "PASS1C-06")
})

test_that("validate_phase error message is context-appropriate", {
  # HUMAN-related path should mention metadata_phase.txt
  expect_error(
    validate_phase("HUMAN-MAIN-TR01/T11/OXYLIPNEG/BATCH1_20240524/PROCESSED_20240524"),
    "metadata_phase.txt"
  )
  # Non-HUMAN path should NOT mention metadata_phase.txt
  err <- tryCatch(
    validate_phase("INVALID/T31/IONPNEG/BATCH1_20190909/PROCESSED_20200205"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "Expected phases")
  expect_false(grepl("metadata_phase.txt", err))
})

test_that("check_viallabel_dmaqc rejects invalid HUMAN-MAIN formats", {
  # Underscore instead of hyphen: HUMAN-MAIN_TR01 (1 hyphen, falls into PRECOVID branch)
  expect_error(
    check_viallabel_dmaqc(
      vl_submitted = "test_vial",
      dmaqc_shipping_info = tempfile(),
      tissue_code = "T11",
      cas = "emory",
      phase = "HUMAN-MAIN_TR01",
      failed_samples = NULL,
      outfile_missed_viallabels = "test"
    ),
    "not a valid HUMAN phase"
  )
  # 2-hyphen but wrong format (e.g., HUMAN-MAIN-XX01)
  expect_error(
    check_viallabel_dmaqc(
      vl_submitted = "test_vial",
      dmaqc_shipping_info = tempfile(),
      tissue_code = "T11",
      cas = "emory",
      phase = "HUMAN-MAIN-XX01",
      failed_samples = NULL,
      outfile_missed_viallabels = "test"
    ),
    "expected format is HUMAN-MAIN-TR##"
  )
})

test_that("set_phase rejects lowercase phase in metadata_phase.txt", {
  # Create a temp directory structure with a lowercase phase in metadata_phase.txt
  tmpdir <- tempfile()
  batch_dir <- file.path(tmpdir, "HUMAN", "T11", "OXYLIPNEG", "BATCH1_20240524")
  proc_dir <- file.path(batch_dir, "PROCESSED_20240524")
  dir.create(proc_dir, recursive = TRUE)
  writeLines("human-main-tr01", file.path(batch_dir, "metadata_phase.txt"))
  expect_error(
    set_phase(input_results_folder = proc_dir, dmaqc_phase2validate = FALSE),
    "must be UPPER CASE"
  )
  unlink(tmpdir, recursive = TRUE)
})

