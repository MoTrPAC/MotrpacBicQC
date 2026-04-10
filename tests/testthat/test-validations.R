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
  # Windows-style backslash paths should also be rejected
  expect_error(validate_phase("D:\\Snow\\HUMAN-MAIN-TR01\\T02\\OXYLIPNEG\\BATCH2_20240808\\PROCESSED_20240808"))
  expect_error(validate_phase("D:\\Snow\\HUMAN-MAIN\\T02\\OXYLIPNEG\\BATCH2_20240808\\PROCESSED_20240808"))
  expect_error(validate_phase("D:\\Snow\\HUMAN_MAIN_TR01\\T02\\OXYLIPNEG\\BATCH2_20240808\\PROCESSED_20240808"))
})

test_that("validate_phase accepts valid folder phases", {
  expect_equal(validate_phase("HUMAN/T11/OXYLIPNEG/BATCH1_20240524/PROCESSED_20240524"), "HUMAN")
  expect_equal(validate_phase("HUMAN-PRECOVID/T10/IONPNEG/BATCH1_20190909/PROCESSED_20200205"), "HUMAN-PRECOVID")
  expect_equal(validate_phase("PASS1A-06/T31/IONPNEG/BATCH1_20190909/PROCESSED_20200205"), "PASS1A-06")
  # Windows-style backslash paths
  expect_equal(validate_phase("D:\\Snow\\HUMAN\\T02\\OXYLIPNEG\\BATCH2_20240808\\PROCESSED_20240808"), "HUMAN")
  expect_equal(validate_phase("D:\\Data\\PASS1A-06\\T31\\IONPNEG\\BATCH1_20190909\\PROCESSED_20200205"), "PASS1A-06")
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
  # Helper: call check_viallabel_dmaqc with a phase and capture the result
  call_vl_dmaqc <- function(phase) {
    check_viallabel_dmaqc(
      vl_submitted = "test_vial",
      dmaqc_shipping_info = tempfile(),
      tissue_code = "T11",
      cas = "emory",
      phase = phase,
      failed_samples = NULL,
      outfile_missed_viallabels = "test"
    )
  }

  # --- Invalid: underscore instead of hyphen (1 hyphen → PRECOVID branch) ---
  expect_error(call_vl_dmaqc("HUMAN-MAIN_TR01"), "not a valid HUMAN phase")
  expect_error(call_vl_dmaqc("HUMAN-MAIN_TR02"), "not a valid HUMAN phase")

  # --- Invalid: 2-hyphen but wrong format ---
  expect_error(call_vl_dmaqc("HUMAN-MAIN-XX01"), "expected format is HUMAN-MAIN-TR##")
  expect_error(call_vl_dmaqc("HUMAN-MAIN-TR1"),  "expected format is HUMAN-MAIN-TR##")   # single digit
  expect_error(call_vl_dmaqc("HUMAN-MAIN-TR001"), "expected format is HUMAN-MAIN-TR##")  # three digits
  expect_error(call_vl_dmaqc("HUMAN-MAIN-tr01"), "expected format is HUMAN-MAIN-TR##")   # lowercase tr
  expect_error(call_vl_dmaqc("HUMAN-FOO-TR01"),  "expected format is HUMAN-MAIN-TR##")   # wrong middle
  expect_error(call_vl_dmaqc("HUMAN-MAIN-"),     "expected format is HUMAN-MAIN-TR##")   # trailing hyphen

  # --- Invalid: 1-hyphen but not HUMAN-PRECOVID ---
  expect_error(call_vl_dmaqc("HUMAN-MAIN"),      "not a valid HUMAN phase")
  expect_error(call_vl_dmaqc("HUMAN-OTHER"),     "not a valid HUMAN phase")
  expect_error(call_vl_dmaqc("HUMAN-precovid"),  "not a valid HUMAN phase")  # lowercase

  # --- Invalid: 3+ hyphens → warning (then errors on file read, but we capture the warning) ---
  expect_warning(
    tryCatch(call_vl_dmaqc("HUMAN-MAIN-TR-01"), error = function(e) NULL),
    "doesn't contain a valid format"
  )

  # --- Invalid: lowercase phase (won't match grepl("HUMAN", eph)) ---
  # lowercase "human" doesn't enter the HUMAN branch at all;
  # it also doesn't match PASS, so it falls through to read the dmaqc file
  # which is an empty tempfile → error from read.delim
  expect_error(suppressWarnings(call_vl_dmaqc("human-main-tr01")))
})

test_that("check_viallabel_dmaqc accepts valid HUMAN phases (parsing only)", {
  # Create a minimal mock dmaqc shipping file so parsing succeeds
  # (will return 0 labels but not error on phase parsing)
  mock_dmaqc <- tempfile(fileext = ".txt")
  writeLines(
    c("vial_label\tbic_tissue_code\tsite_code\tphase\ttranche\tanimal_age",
      "VL001\tT11\temory\tHUMAN\t00\tNA"),
    mock_dmaqc
  )

  call_vl_dmaqc_valid <- function(phase) {
    check_viallabel_dmaqc(
      vl_submitted = "VL001",
      dmaqc_shipping_info = mock_dmaqc,
      tissue_code = "T11",
      cas = "emory",
      phase = phase,
      failed_samples = NULL,
      outfile_missed_viallabels = "test",
      return_n_issues = TRUE,
      verbose = FALSE
    )
  }

  # Valid: plain HUMAN (0 hyphens)
  expect_no_error(call_vl_dmaqc_valid("HUMAN"))

  # Valid: HUMAN-PRECOVID (1 hyphen)
  expect_no_error(call_vl_dmaqc_valid("HUMAN-PRECOVID"))

  # Valid: HUMAN-MAIN-TR## (2 hyphens) — various tranche numbers
  expect_no_error(call_vl_dmaqc_valid("HUMAN-MAIN-TR01"))
  expect_no_error(call_vl_dmaqc_valid("HUMAN-MAIN-TR02"))
  expect_no_error(call_vl_dmaqc_valid("HUMAN-MAIN-TR10"))
  expect_no_error(call_vl_dmaqc_valid("HUMAN-MAIN-TR99"))

  unlink(mock_dmaqc)
})

test_that("set_phase validates metadata_phase.txt content", {
  # Helper: create temp folder structure and write metadata_phase.txt
  setup_phase_dir <- function(phase_content) {
    tmpdir <- tempfile()
    batch_dir <- file.path(tmpdir, "HUMAN", "T11", "OXYLIPNEG", "BATCH1_20240524")
    proc_dir <- file.path(batch_dir, "PROCESSED_20240524")
    dir.create(proc_dir, recursive = TRUE)
    writeLines(phase_content, file.path(batch_dir, "metadata_phase.txt"))
    list(proc_dir = proc_dir, tmpdir = tmpdir)
  }

  call_set_phase <- function(phase_content) {
    dirs <- setup_phase_dir(phase_content)
    on.exit(unlink(dirs$tmpdir, recursive = TRUE))
    set_phase(input_results_folder = dirs$proc_dir, dmaqc_phase2validate = FALSE, verbose = FALSE)
  }

  # --- Reject lowercase ---
  expect_error(call_set_phase("human-main-tr01"), "must be UPPER CASE")
  expect_error(call_set_phase("Human-Main-TR01"), "must be UPPER CASE")
  expect_error(call_set_phase("pass1a-06"),       "must be UPPER CASE")

  # --- Reject underscore variant ---
  expect_error(call_set_phase("HUMAN-MAIN_TR01"), "Expected format: HUMAN-MAIN-TR##")

  # --- Reject incomplete HUMAN-MAIN ---
  expect_error(call_set_phase("HUMAN-MAIN"),      "Expected format: HUMAN-MAIN-TR##")
  expect_error(call_set_phase("HUMAN-MAIN-TR"),   "Expected format: HUMAN-MAIN-TR##")
  expect_error(call_set_phase("HUMAN-MAIN-TR1"),  "Expected format: HUMAN-MAIN-TR##")
  expect_error(call_set_phase("HUMAN-MAIN-TR001"), "Expected format: HUMAN-MAIN-TR##")

  # --- Accept valid HUMAN-MAIN-TR## ---
  expect_equal(call_set_phase("HUMAN-MAIN-TR01"), "HUMAN-MAIN-TR01")
  expect_equal(call_set_phase("HUMAN-MAIN-TR02"), "HUMAN-MAIN-TR02")
  expect_equal(call_set_phase("HUMAN-MAIN-TR10"), "HUMAN-MAIN-TR10")

  # --- Accept valid PASS phases ---
  dirs <- setup_phase_dir("PASS1A-06")
  # This will fail because the folder is HUMAN/ but phase says PASS — that's fine,
  # we're testing that the uppercase and format check passes
  expect_equal(
    set_phase(input_results_folder = dirs$proc_dir, dmaqc_phase2validate = FALSE, verbose = FALSE),
    "PASS1A-06"
  )
  unlink(dirs$tmpdir, recursive = TRUE)
})

test_that("validate_assay works with Windows backslash paths", {
  expect_equal(validate_assay("D:\\Snow\\HUMAN\\T02\\OXYLIPNEG\\BATCH2_20240808\\PROCESSED_20240808"), "OXYLIPNEG")
  expect_equal(validate_assay("D:\\Data\\PASS1A-06\\T31\\IONPNEG\\BATCH1_20190909\\PROCESSED_20200205"), "IONPNEG")
  expect_equal(validate_assay("C:\\Users\\lab\\HUMAN\\T10\\PROT_AC\\BATCH1_20230620\\PROCESSED_20230620"), "PROT_AC")
  expect_equal(validate_assay("D:\\PASS1B-06\\T58\\LAB_CK\\BATCH1_20230620\\RESULTS_20230620"), "LAB_CK")
})

test_that("validate_batch works with Windows backslash paths", {
  expect_equal(
    validate_batch("D:\\Snow\\HUMAN\\T02\\OXYLIPNEG\\BATCH2_20240808\\PROCESSED_20240808"),
    "D:\\Snow\\HUMAN\\T02\\OXYLIPNEG\\BATCH2_20240808\\"
  )
  expect_equal(
    validate_batch("D:\\Data\\PASS1A-06\\T31\\IONPNEG\\BATCH1_20190909\\PROCESSED_20200205"),
    "D:\\Data\\PASS1A-06\\T31\\IONPNEG\\BATCH1_20190909\\"
  )
  # PROT_AC legacy on Windows
  expect_equal(
    validate_batch("C:\\broad\\PASS1A-06\\T58\\PROT_AC\\BATCH_20190828\\PROCESSED_20200101"),
    "C:\\broad\\PASS1A-06\\T58\\PROT_AC\\BATCH_20190828\\"
  )
})

test_that("get_full_path2batch works with Windows backslash paths", {
  expect_equal(
    get_full_path2batch("D:\\Snow\\HUMAN\\T02\\OXYLIPNEG\\BATCH2_20240808\\PROCESSED_20240808"),
    "D:\\Snow\\HUMAN\\T02\\OXYLIPNEG\\BATCH2_20240808\\"
  )
  expect_equal(
    get_full_path2batch("D:\\Data\\HUMAN\\T10\\LAB_CK\\BATCH1_20230620\\RESULTS_20230620"),
    "D:\\Data\\HUMAN\\T10\\LAB_CK\\BATCH1_20230620\\"
  )
})

