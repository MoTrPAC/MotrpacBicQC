
#' @title Check analyte metadata file
#'
#' @description Checks whether the metadata_analyte file follows the required guidelines.
#' @param df (data.frame) The metadata_analyte data frame to check.
#' @param return_n_issues (logical) If `TRUE`, returns the number of issues found.
#' @param validate_uniprot (logical) If `TRUE`, checks if all Uniprot IDs are valid by connecting to the Uniprot database. Note: This may take several minutes depending on the number of IDs.
#' @param verbose (logical) If `TRUE` (default), displays messages during the checking process.
#' @return (int) The number of issues identified if `return_n_issues` is `TRUE`.
#' @examples
#' \dontrun{
#' check_metadata_analyte(df = metadata_analyte)
#' }
#' @export
check_metadata_analyte <- function(df,
                                   return_n_issues = FALSE,
                                   validate_uniprot = FALSE,
                                   verbose = TRUE) {
  
  # Initialize issue count
  ic <- 0
  
  # Ensure required columns are present
  df <- filter_required_columns(df = df,
                                type = "labanalytes",
                                verbose = verbose)
  
  # Check 'analyte_name' column
  if ("analyte_name" %in% colnames(df)) {
    # Check for uniqueness
    if (length(unique(df$analyte_name)) != nrow(df)) {
      duplis_details <- df$analyte_name[duplicated(df$analyte_name)]
      duplis <- length(unique(duplis_details))
      if (verbose) message("   - (-) `analyte_name` non-unique values detected: ", duplis)
      if (verbose) message("\t\t - ", paste(unique(duplis_details), collapse = "\n\t\t - "))
      ic <- ic + 1
    } else {
      if (verbose) message("  + (+) `analyte_name` unique values: OK")
    }
    # Check for missing values
    if (check_missing_values(df, "analyte_name")) {
      if (verbose) message("   - (-) `analyte_name`: NA values detected: FAIL")
      ic <- ic + 1
    }
  } else {
    if (verbose) message("   - (-) `analyte_name` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Check 'uniprot_entry' column
  if ("analyte_id" %in% colnames(df)) {
    # Check for uniqueness (duplicates may be acceptable)
    if (length(unique(df$analyte_id)) != nrow(df)) {
      duplis_details <- df$analyte_id[duplicated(df$analyte_id)]
      duplis <- length(unique(duplis_details))
      if (verbose) message("   - (-) `analyte_id` non-unique values detected (n duplications = ", duplis, "). FAIL")
      if (verbose) message("\t\t - ", paste(unique(duplis_details), collapse = "\n\t\t - "))
    } else {
      if (verbose) message("  + (+) `analyte_id` unique values: OK")
    }
    # Validate Uniprot IDs if requested
    # if (validate_uniprot) {
    #   if (verbose) message("  + Validating `uniprot_entry` IDs with the Uniprot database. Please wait...")
    #   all_valid <- validate_uniprot_ids_with_uniprot(ids = unique(df$uniprot_entry))
    #   if (all_valid) {
    #     if (verbose) message("  + (+) All `uniprot_entry` IDs are valid: OK")
    #   } else {
    #     if (verbose) message("  + (-) Some `uniprot_entry` IDs are invalid: FAIL")
    #     ic <- ic + 1
    #   }
    # }
    # # Check for missing values
    # if (check_missing_values(df, "uniprot_entry")) {
    #   if (verbose) message("   - (-) `uniprot_entry`: NA values detected: FAIL")
    #   ic <- ic + 1
    # }
  } else {
    if (verbose) message("   - (-) `analyte_id` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Check 'database_id' column
  if ("database_id" %in% colnames(df)) {
    if (check_missing_values(df, "database_id")) {
      if (verbose) message("   - (-) `database_id`: NA values detected: FAIL")
      ic <- ic + 1
    }else{
      if (verbose) message("  + (+) `database_id` values are valid: OK")
    }
  } else {
    if (verbose) message("   - (-) `assay_name` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Check 'assay_name' column
  if ("assay_name" %in% colnames(df)) {
    # Check for uniqueness (duplicates may be acceptable)
    if (length(unique(df$assay_name)) != nrow(df)) {
      duplis_details <- df$assay_name[duplicated(df$assay_name)]
      duplis <- length(unique(duplis_details))
      if (verbose) message("   - ( ) `assay_name` non-unique values detected (n duplications = ", duplis, "). This is acceptable.")
      if (verbose) message("\t\t - ", paste(unique(duplis_details), collapse = "\n\t\t - "))
    } else {
      if (verbose) message("  + (+) `assay_name` unique values: OK")
    }
    # Check for missing values
    if (check_missing_values(df, "assay_name")) {
      if (verbose) message("   - (-) `assay_name`: NA values detected: FAIL")
      ic <- ic + 1
    }
  } else {
    if (verbose) message("   - (-) `assay_name` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Check 'unit' column
  if ("unit" %in% colnames(df)) {
    # Check for uniqueness (duplicates may be acceptable)
    if (length(unique(df$unit)) != nrow(df)) {
      duplis_details <- df$unit[duplicated(df$unit)]
      duplis <- length(unique(duplis_details))
      if (verbose) message("   - ( ) `unit` non-unique values detected (n duplications = ", duplis, "). This is acceptable.")
      if (verbose) message("\t\t - ", paste(unique(duplis_details), collapse = "\n\t\t - "))
    } else {
      if (verbose) message("  + (+) `unit` unique values: OK")
    }
    # Check for missing values
    if (check_missing_values(df, "unit")) {
      if (verbose) message("   - (-) `unit`: NA values detected: FAIL")
      ic <- ic + 1
    }
  } else {
    if (verbose) message("   - (-) `unit` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Return the number of issues if requested
  if (return_n_issues) return(ic)
  
} # End of check_metadata_analyte


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Check LAB metadata samples file
#'
#' @description Checks whether the metadata_sample file for LAB assays follows the required guidelines.
#' @param df (data.frame) LAB metadata samples data
#' @param return_n_issues (logical) If `TRUE`, returns the number of issues identified.
#' @param verbose (logical) If `TRUE` (default), displays messages during the checking process.
#' @return (int) Number of issues identified if `return_n_issues` is `TRUE`.
#' @examples
#' \dontrun{
#' check_metadata_samples_lab(df = metadata_sample_named_CK_plasma)
#' }
#' @export
check_metadata_samples_lab <- function(df,
                                       return_n_issues = FALSE,
                                       verbose = TRUE) {
  
  # Initialize issue count
  ic <- 0
  
  # Filter only expected columns
  df <- filter_required_columns(df = df,
                                type = "labsamples",
                                verbose = verbose)
  
  # Check 'sample_id' column
  if ("sample_id" %in% colnames(df)) {
    if (length(unique(df$sample_id)) != nrow(df)) {
      if (verbose) message("   - (-) `sample_id`: Non-unique values detected: FAIL")
      ic <- ic + 1
    } else {
      if (verbose) message("  + (+) `sample_id` unique values: OK")
    }
  } else {
    if (verbose) message("   - (-) `sample_id` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Define expected sample types
  expected_sample_types <- c("Sample", "QC-Pooled", "QC-Reference", "QC-Blank",
                             "QC-Identification", "QC-InternalStandard", "QC-PreRun",
                             "QC-ExternalStandard", "QC-DriftCorrection", "QC-ReCAS", 
                             "QC-PlateControl")
  
  # Check 'sample_type' column
  if ("sample_type" %in% colnames(df)) {
    undefined_types <- setdiff(df$sample_type, expected_sample_types)
    if (length(undefined_types) > 0) {
      if (verbose) message("   - (-) Undefined `sample_type` values detected: FAIL")
      if (verbose) message("\t\t - ", paste(undefined_types, collapse = "\n\t\t - "))
      ic <- ic + 1
    } else {
      if (verbose) message("  + (+) `sample_type` values are valid: OK")
    }
  } else {
    if (verbose) message("   - (-) `sample_type` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Check 'sample_order' column
  if ("sample_order" %in% colnames(df)) {
    if (!all(is.numeric(df$sample_order))) {
      if (verbose) message("   - (-) `sample_order` contains non-numeric values: FAIL")
      non_numeric_values <- df$sample_order[!is.numeric(df$sample_order)]
      if (verbose) message("\t\t - ", paste(non_numeric_values, collapse = "\n\t\t - "))
      ic <- ic + 1
    } else {
      if (verbose) message("  + (+) `sample_order` is numeric: OK")
    }
  } else {
    if (verbose) message("   - (-) `sample_order` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Check 'raw_file' column
  if ("raw_file" %in% colnames(df)) {
    if (any(is.na(df$raw_file) | df$raw_file == "")) {
      if (verbose) message("   - (-) `raw_file` contains missing or empty values: FAIL")
      ic <- ic + 1
    } else {
      if (verbose) message("  + (+) `raw_file` values are valid: OK")
    }
  } else {
    if (verbose) message("   - (-) `raw_file` column missing: FAIL")
    ic <- ic + 1
  }
  
  if("extraction_date" %in% colnames(df)){
    if(any(is.na(df$extraction_date))){
      if(verbose) message("   - (-) `extraction_date` has NA values: FAIL")
      ic <- ic + 1
    }else{
      icdate <- validate_yyyymmdd_dates(df = df, date_column = "extraction_date", verbose = verbose)
      ic <- ic + icdate
    }
  }else{
    if(verbose) message("   - (-) `extraction_date` column missed: FAIL")
    ic <- ic + 1
  }
  
  if("acquisition_date" %in% colnames(df)){
    if(any(is.na(df$acquisition_date))){
      if(verbose) message("   - (-) `acquisition_date` has NA values: FAIL")
      ic <- ic + 1
    }else{
      if( any(grepl(":", df$acquisition_date)) ){
        if(verbose) message("  + (i) Assuming `acquisition_date` is in `MM/DD/YYYY HH:MM:SS AM/PM` format. Validating:")
        icdt <- validate_dates_times(df = df, column_name = "acquisition_date", verbose = verbose)
      }else{
        icdate <- validate_yyyymmdd_dates(df = df, date_column = "acquisition_date", verbose = verbose)
        ic <- ic + icdate
      }
    }
  }else{
    if(verbose) message("   - (-) `acquisition_date` column missed: FAIL")
    ic <- ic + 1
  }
  
  # Return the number of issues if requested
  if (return_n_issues) return(ic)
  
} # End of check_metadata_samples_lab

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Load and Process LAB Batch Data
#'
#' @description
#' This function loads LAB batch data from the specified input directory.
#' It performs quality checks on the data and loads specific files related to
#' LAB data, including metadata for analytes and samples, and the results file.
#' It also integrates validation checks and warns if there are too many issues
#' identified in the data.
#'
#' @param input_results_folder A string representing the path to the folder
#' containing LAB batch data to be loaded and processed.
#' @param verbose Logical; if `TRUE`, prints detailed messages
#' during the loading process.
#'
#' @return A list containing data frames for metadata of analytes (`m_a`),
#' metadata of samples (`m_s`), and LAB results (`r_o`). If certain files
#' are not available, the corresponding entries in the list will be `NULL`.
#' @examples
#' \dontrun{
#' list_of_df <- load_lab_batch(input_results_folder = "/path/to/PROCESSED_YYYYMMDD/")
#' }
#' @export
load_lab_batch <- function(input_results_folder,
                           verbose = TRUE) {
  
  m_a <- m_s <- r_o <- NULL
  
  # Validations ----
  phase <- validate_phase(input_results_folder)
  processfolder <- validate_processFolder(input_results_folder)
  assay <- validate_assay(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)
  
  # Perform validation using validate_lab
  total_issues <- validate_lab(input_results_folder = input_results_folder,
                               cas = "duke",
                               return_n_issues = TRUE,
                               verbose = FALSE)
  
  if (total_issues > 0) {
    message("\n\tWARNING!!! Too many issues identified (", total_issues, ").\n",
            "\tThis batch should not be processed until the issues are resolved.")
  }
  
  vial_label <- NA
  qc_samples <- NA
  
  # Load LAB data ----
  if (verbose) message("# LOAD LAB BATCH")
  if (verbose) message("+ Site: Duke, Chris Newgard group")
  if (verbose) message("+ Folder: `", paste0(input_results_folder), "`")
  
  # Load metadata_analyte file ----
  if (verbose) message("\n## Loading `metadata_analyte` file\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata_analyte_.*\\.txt$",
                     verbose = verbose)
  f_ma <- lista$flag
  if (f_ma) {
    m_a_f <- lista$filename
    m_a <- filter_required_columns(df = lista$df,
                                  type = "labanalytes",
                                  verbose = verbose)
  } else {
    if (verbose) message("   - (-) `metadata_analyte` file not available")
  }
  
  # Load metadata_sample file ----
  if (verbose) message("\n## Loading `metadata_sample` file\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata_sample_.*\\.txt$",
                     verbose = verbose)
  f_ms <- lista$flag
  if (f_ms) {
    m_s_f <- lista$filename
    m_s <- filter_required_columns(df = lista$df,
                                   type = "labsamples",
                                   verbose = verbose)
  } else {
    if (verbose) message("   - (-) `metadata_sample` file not available")
  }
  
  # Load results file ----
  if (verbose) message("\n## Loading `results` file\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results_.*\\.txt$",
                     verbose = verbose)
  f_r <- lista$flag
  if (f_r) {
    r_f <- lista$filename
    r_o <- lista$df
  } else {
    if (verbose) message("   - (-) `results` file not available")
  }
  
  # RETURN list of data frames ----
  listdf <- list("m_a" = m_a,
                 "m_s" = m_s,
                 "r_o" = r_o)
  
  return(listdf)
}


#' Validate LAB Assay Data
#'
#' This function validates LAB assay data based on a set of criteria, including
#' folder structure, metadata, and file contents. It supports optional
#' functionalities like creating a PDF report and DMAQC validation
#' (data only available at the BIC).
#'
#' @param input_results_folder A string representing the path to the folder
#'        containing LAB assay results to be validated.
#' @param cas A character string indicating the CAS number.
#' @param return_n_issues Logical; if `TRUE`, the function returns the number of
#'        detected issues.
#' @param full_report Logical; if `TRUE`, generates a full report of the
#'        validation process.
#' @param f_proof Logical; if `TRUE`, generates proof plots for data validation.
#' @param printPDF Logical; if `TRUE` and `f_proof` is `TRUE`, saves the plots
#'        to a PDF file. If so, provide the desired path to output
#'        the PDF file in the argument `out_qc_folder`.
#' @param out_qc_folder Optional; a string specifying the path to the folder
#'        where output PDF should be saved (only if `printPDF = TRUE`).
#'        Default: current working directory.
#' @param dmaqc_shipping_info (character) File path to the DMAQC file.
#'        Only the BIC can use this argument.
#' @param dmaqc_phase2validate (character) Provide phase to validate. This argument
#'        is not required since it should be extracted from the input folder or from the
#'        new required file `metadata_phase.txt`. Please ignore. However, if this argument is provided,
#'        it will take priority and this will be the phase.
#' @param validate_uniprot Logical; if `TRUE`, validates against the UniProt
#'        database.
#' @param verbose Logical; if `TRUE`, prints detailed messages during validation.
#'
#' @return Depending on the settings, this function may return the number of
#' issues found, generate reports or plots, or simply perform the
#' validation without returning anything.
#' @examples
#' \dontrun{
#' validate_lab("/path/to/results", cas = "broad_rg", return_n_issues = TRUE)
#' }
#' @export
validate_lab <- function(input_results_folder,
                         cas,
                         return_n_issues = FALSE,
                         full_report = FALSE,
                         f_proof = FALSE,
                         printPDF = FALSE,
                         out_qc_folder = NULL,
                         dmaqc_shipping_info = NULL,
                         dmaqc_phase2validate = FALSE,
                         validate_uniprot = FALSE,
                         verbose = TRUE) {
  
  analyte_name = plot_basic_lab_qc = sample_id = sample_order = value = NULL
  
  # Initialize issue counts
  ic <- 0       # Total issues
  ic_m_a <- 0   # Issues in metadata analyte file
  ic_m_s <- 0   # Issues in metadata sample file
  ic_r <- 0     # Issues in results file
  ic_man <- 0   # Issues in manifest file
  
  # Placeholders for vial label and QC samples counts
  vial_label <- NA
  qc_samples <- NA
  
  # Placeholder for DMAQC validation issues
  if (is.null(dmaqc_shipping_info)) {
    ic_vl <- "missed"
  } else {
    ic_vl <- NA
  }
  
  # Validate folder structure and extract metadata
  validate_cas(cas = cas)
  processfolder <- validate_processFolder(input_results_folder)
  assay <- validate_assay(input_results_folder)
  phase <- validate_phase(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)
  batch_folder <- validate_batch(input_results_folder)
  
  input_results_folder <- normalizePath(input_results_folder)
  
  # Extract short folder path for reporting
  input_folder_short <- regmatches(input_results_folder, regexpr("(HUMAN|PASS).*RESULTS_[0-9]{8}", input_results_folder))
  if (purrr::is_empty(input_folder_short)) {
    if (verbose) message("\nThe RESULTS_YYYYMMDD folder full path is not correct. Example:")
    if (verbose) message("/full/path/to/folder/HUMAN/T02/LAB_CK/BATCH1_20230620/PROCESSED_20230620")
    stop("Input folder not according to guidelines")
  }
  
  if (verbose) message("# LAB Assay QC Report\n\n")
  if (verbose) message("+ Site: ", cas)
  if (verbose) message("+ Folder: `", paste0(input_folder_short), "`")
  
  # Check for metadata phase file----
  is_mp <- check_metadata_phase_file(input_results_folder = input_results_folder, verbose = verbose)
  if (!is_mp) {
    ic <- ic + 1
  }
  
  # Set phase----
  dmaqc_phase2validate <- set_phase(input_results_folder = input_results_folder,
                                    dmaqc_phase2validate = dmaqc_phase2validate,
                                    verbose = verbose)
  
  phase2file <- generate_phase_details(phase_metadata = dmaqc_phase2validate)
  
  # QC metadata_analyte file----
  if (verbose) message("\n## QC `metadata_analyte` file\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata_analyte_.*\\.txt$",
                     verbose = verbose)
  f_ma <- lista$flag
  if (f_ma) {
    m_a_f <- lista$filename
    m_a <- lista$df
    ic_m_a <- check_metadata_analyte(df = m_a,
                                     return_n_issues = TRUE,
                                     validate_uniprot = validate_uniprot,
                                     verbose = verbose)
  } else {
    if (verbose) message("   - (-) `metadata_analyte` file not available: FAIL")
    ic_m_a <- 10
  }
  
  # QC metadata_sample file----
  if (verbose) message("\n## QC `metadata_sample` file\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata_sample_.*\\.txt$",
                     verbose = verbose)
  f_ms <- lista$flag
  if (f_ms) {
    m_s_f <- lista$filename
    m_s <- lista$df
    ic_m_s <- check_metadata_samples_lab(df = m_s,
                                         return_n_issues = TRUE,
                                         verbose = verbose)
    # Extract the number of samples
    if (!is.null(m_s)) {
      if (all(c("sample_id", "sample_type") %in% colnames(m_s))) {
        vial_label <- sum(m_s$sample_type == "Sample")
        qc_samples <- sum(m_s$sample_type != "Sample")
      }
    }
  } else {
    if (verbose) message("   - (-) `metadata_sample` file not available: FAIL")
    ic_m_s <- 10
  }
  
  # QC results file-----
  if (verbose) message("\n## QC `results` file\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results_.*\\.txt$",
                     verbose = verbose)
  f_r <- lista$flag
  if (f_r) {
    r_f <- lista$filename
    r_o <- lista$df
    ic_r <- check_results_assays(df = r_o,
                          return_n_issues = TRUE,
                          assay_type = "lab",
                          verbose = verbose)
  } else {
    if (verbose) message("   - (-) `results` file not available: FAIL")
    ic_r <- 10
  }
  
  # Cross-file validation------
  if (verbose) message("\n## Cross-File Validation\n")
  if (ic_m_a == 0 && ic_m_s == 0 && ic_r == 0) {
    ic_c_f_v <- check_crossfile_validation(r_o = r_o,
                                           m_s = m_s,
                                           m_p = m_a,
                                           return_n_issues = TRUE,
                                           assay_type = "lab",
                                           verbose = verbose)
  } else {
    if (verbose) message("   - (-) Cross-file validation is not possible due to previous errors.")
    ic <- ic + 1
  }
  
  # QC Plots-----
  if (f_proof) {
    if (verbose) message("\n## QC Plots\n")
    output_prefix <- paste0(cas, ".", tolower(phase2file), ".", tissue_code, ".", tolower(assay), ".", tolower(processfolder))
    output_prefix <- gsub("\\|", "_", output_prefix)
    
    if (f_r & f_ms & f_ma) {
      if (nrow(r_o) != 0) {
        # Prepare data for plotting
        results_long <- r_o %>%
          tidyr::pivot_longer(cols = -analyte_name,
                              names_to = "sample_id",
                              values_to = "value")
        results_long <- merge(m_s, results_long, by = "sample_id")
        results_long <- results_long %>%
          filter(!is.na(value), value != 0) %>%
          arrange(sample_order) %>%
          mutate(sample_id_ordered = factor(sample_id, levels = unique(sample_id)))
        
        # Merge analyte metadata
        r_p <- merge(m_a, r_o, by = "analyte_name")
        
        # Create output folder-------
        if (is.null(out_qc_folder)){
          out_qc_folder <- getwd()
        }else{
          out_qc_folder <- create_folder(out_qc_folder)
        }
        
        # Generate plots
        plot_basic_lab_qc(results = r_p,
                          results_long = results_long,
                          out_qc_folder = out_qc_folder,
                          output_prefix = output_prefix,
                          printPDF = printPDF,
                          verbose = verbose)
      } else {
        if (verbose) message("  (-) QC plots are not possible: not enough analytes.")
      }
    } else {
      if (verbose) message("  (-) QC plots are not possible: critical datasets are missing.")
    }
  }
  
  # Manifest file check-----
  if (verbose) message("\n## QC `file_manifest_YYYYMMDD.csv` (required)\n")
  
  file_manifest <- list.files(normalizePath(batch_folder),
                              pattern = "file_manifest.*csv",
                              ignore.case = TRUE,
                              full.names = TRUE,
                              recursive = TRUE)
  
  if (length(file_manifest) == 0) {
    f_man <- FALSE
    if (verbose) message("   - (-) `file_manifest_YYYYMMDD.csv` file not available: FAIL")
    ic_man <- ic_man + 6
    ic <- ic + 1
  } else {
    first_manifest <- sort(basename(file_manifest), decreasing = TRUE)[1]
    file_manifest <- file_manifest[basename(file_manifest) == first_manifest]
    f_man <- TRUE
  }
  
  if (f_man) {
    manifest <- read.csv(file_manifest)
    required_columns <- c("file_name", "md5")
    if (all(required_columns %in% colnames(manifest))) {
      if (verbose) message("  + (+) `file_name`, `md5` columns available in manifest file.")
      
      # Normalize file paths
      manifest$file_name <- gsub("\\\\", "/", manifest$file_name)
      manifest$file_base <- basename(manifest$file_name)
      
      # Check for metadata_analyte file
      if (f_ma) {
        analyte_file_in_manifest <- basename(manifest$file_name[grepl("metadata_analyte", manifest$file_name)])
        to_check <- basename(m_a_f)
        if (to_check %in% analyte_file_in_manifest) {
          if (verbose) message("  + (+) `metadata_analyte` file included in manifest: OK")
        } else {
          if (verbose) message("   - (-) `metadata_analyte` file not included in manifest: FAIL")
          ic_man <- ic_man + 1
        }
      }
      
      # Check for metadata_sample file
      if (f_ms) {
        sample_file_in_manifest <- basename(manifest$file_name[grepl("metadata_sample", manifest$file_name)])
        to_check <- basename(m_s_f)
        if (to_check %in% sample_file_in_manifest) {
          if (verbose) message("  + (+) `metadata_sample` file included in manifest: OK")
        } else {
          if (verbose) message("   - (-) `metadata_sample` file not included in manifest: FAIL")
          ic_man <- ic_man + 1
        }
      }
      
      # Check for results file
      if (f_r) {
        results_file_in_manifest <- basename(manifest$file_name[grepl("results", manifest$file_name)])
        to_check <- basename(r_f)
        if (to_check %in% results_file_in_manifest) {
          if (verbose) message("  + (+) `results` file included in manifest: OK")
        } else {
          if (verbose) message("   - (-) `results` file not included in manifest: FAIL")
          ic_man <- ic_man + 1
        }
      }
      
      # Check for missing MD5 values
      if (any(is.na(manifest$md5))) {
        if (verbose) message("   - (-) `md5` column contains NA values: FAIL")
        ic_man <- ic_man + 1
      }
    } else {
      if (verbose) message("   - (-) Not all required columns are available in manifest: FAIL")
      ic_man <- ic_man + 1
    }
  }
  
  # DMAQC validation-----
  if (verbose) message("\n## DMAQC Validation\n")
  
  failed_samples <- check_failedsamples(input_results_folder = batch_folder, verbose = verbose)
  
  # Validate vial labels from DMAQC
  if (is.na(ic_vl)) {
    if (f_ms) {
      vl_results <- m_s$sample_id[m_s$sample_type == "Sample"]
      outfile_missed_viallabels <- paste0(cas, ".", tolower(phase2file), ".", tissue_code, ".", tolower(assay), ".", tolower(processfolder))
      outfile_missed_viallabels <- gsub("\\|", "_", outfile_missed_viallabels)
      ic_vl <- check_viallabel_dmaqc(vl_submitted = vl_results,
                                     tissue_code = tissue_code,
                                     cas = cas,
                                     phase = dmaqc_phase2validate,
                                     failed_samples = failed_samples,
                                     dmaqc_shipping_info = dmaqc_shipping_info,
                                     out_qc_folder = out_qc_folder,
                                     outfile_missed_viallabels = outfile_missed_viallabels,
                                     return_n_issues = TRUE,
                                     verbose = verbose)
    }
  }
  
  # Final reporting
  if (ic > 4) {
    message("\nTOTAL NUMBER OF CRITICAL ERRORS: ", ic, "\n")
    message("WARNING: Too many errors. Please revise the input folder.")
  }
  
  batchversion <- stringr::str_extract(string = input_results_folder, pattern = "BATCH.*_[0-9]+/RESULTS_[0-9]+")
  
  qc_date <- format(Sys.time(), "%Y%m%d_%H%M%S")
  t_name <- bic_animal_tissue_code$bic_tissue_name[bic_animal_tissue_code$bic_tissue_code == tissue_code]
  
  if (return_n_issues) {
    # Only include DMAQC issues in total if dmaqc_shipping_info was provided
    if (!is.null(dmaqc_shipping_info) && is.numeric(ic_vl)) {
      total_issues <- sum(ic, ic_man, ic_m_a, ic_m_s, ic_r, ic_vl, na.rm = TRUE)
    } else {
      total_issues <- sum(ic, ic_man, ic_m_a, ic_m_s, ic_r, na.rm = TRUE)
    }
    
    if (verbose) message("\nTOTAL NUMBER OF ISSUES: ", total_issues, "\n")
    if (full_report) {
      reports <- data.frame(
        cas = cas,
        phase = dmaqc_phase2validate,
        tissue = tissue_code,
        t_name = t_name,
        assay = assay,
        version = batchversion,
        vial_label = vial_label,
        qc_samples = qc_samples,
        dmaqc_valid = ic_vl,
        critical_issues = ic,
        manifest = ic_man,
        m_analyte = ic_m_a,
        m_sample = ic_m_s,
        results = ic_r,
        qc_date = qc_date
      )
      return(reports)
    } else {
      return(total_issues)
    }
  }
} # End of validate_lab


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Write LAB Data Release
#'
#' @description Write out LAB data releases. This function doesn't check whether
#' data has been submitted according to guidelines; it assumes that the data
#' has already been validated and is ready for release.
#'
#' @param input_results_folder (character) Path to the `PROCESSED_YYYYMMDD` folder
#' containing the LAB data.
#' @param folder_name (character) Output folder name. Default is `"motrpac_release"`.
#' @param folder_root (character) Absolute path where the output folder will be created.
#' Default is the current working directory.
#' @param version_file (character) File version number (e.g., `"v1.0"`). Default is `"v1.0"`.
#' @param verbose (logical) If `TRUE` (default), displays messages during processing.
#' @return Creates the BIC release folder and file structure:
#' `PHASE/OMICS/TCODE_NAME/ASSAY/`, and writes out files with names:
#' `motrpac_phase-code_tissuecode_assay_file-details-version.txt`,
#' where `file-details` can be:
#' - `metadata-analyte`, 
#' - `metadata-samples`,
#' - `results`
#' @examples
#' \dontrun{
#' write_lab_releases(
#'    input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD/")
#' }
#' @export
write_lab_releases <- function(input_results_folder,
                               folder_name = "motrpac_release",
                               folder_root = NULL,
                               version_file = "v1.0",
                               verbose = TRUE) {
  
  # Get names from input_results_folder ------
  assay <- validate_assay(input_results_folder)
  phase <- validate_phase(input_results_folder)
  phase_metadata <- set_phase(input_results_folder = input_results_folder, 
                              dmaqc_phase2validate = FALSE, 
                              verbose = FALSE)
  phase_details <- generate_phase_details(phase_metadata)
  tissue_code <- validate_tissue(input_results_folder)
  
  # Get folder_tissue from bic_animal_tissue_code
  folder_tissue <- bic_animal_tissue_code$tissue_name_release[
    which(bic_animal_tissue_code$bic_tissue_code == tissue_code)
  ]
  
  # Get folder_assay from assay_codes
  if (length(assay_codes$assay_code[assay_codes$submission_code == assay]) == 1) {
    folder_assay <- assay_codes$assay_code[assay_codes$submission_code == assay]
  } else {
    stop("ASSAY code ", assay, " not available in `assay_codes`")
  }
  
  if (verbose) message("+ Writing out: ", phase_details, " ", tissue_code, " ", assay, " files", appendLF = FALSE)
  
  # Load LAB datasets ----
  lab_df <- load_lab_batch(input_results_folder = input_results_folder, verbose = FALSE)
  
  # Create output folder -------
  if (is.null(folder_root)) {
    folder_root <- getwd()
  } else {
    folder_root <- create_folder(folder_root)
  }
  
  # Exception for PASS1C-06: the main folder is pass1a
  if (phase_details == "pass1c-06") {
    phase_folder_release <- "pass1a-06"
  } else {
    phase_folder_release <- phase_details
  }
  
  # Decide on the OMICS folder for LAB data
  # LAB assays are typically classified under "clinical-chemistry"
  omics_folder <- "clinical-chemistry"
  
  # Construct the output folder path
  output_folder <- file.path(folder_root, folder_name, phase_folder_release, "results", omics_folder, folder_tissue, folder_assay)
  
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Create shared file name
  file_name_shared <- paste0("motrpac_",
                             phase_details, "_",
                             folder_tissue, "_",
                             folder_assay)
  
  # Create and write FILES -----
  metadata_analyte <- file.path(output_folder, paste0(file_name_shared, "_metadata-analyte_", version_file, ".txt"))
  metadata_samples <- file.path(output_folder, paste0(file_name_shared, "_metadata-samples_", version_file, ".txt"))
  results <- file.path(output_folder, paste0(file_name_shared, "_results_", version_file, ".txt"))
  
  write.table(lab_df$m_a, metadata_analyte, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(lab_df$m_s, metadata_samples, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(lab_df$r, results, row.names = FALSE, sep = "\t", quote = FALSE)
  
  if (verbose) message("...done!")
}

