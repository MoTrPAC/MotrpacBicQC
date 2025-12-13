# VALIDATIONS


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Cross-validate assay files (Olink and LAB)
#'
#' @description Checks that values from the results file are available in both the
#' metadata analyte/protein and metadata sample files.
#' @param r_o (data.frame) Results data frame.
#' @param m_s (data.frame) Metadata sample data frame.
#' @param m_p (data.frame) Metadata analyte/protein data frame.
#' @param assay_type (character) The type of assay, either `"olink"` or `"lab"`.
#' @param return_n_issues (logical) If `TRUE`, returns the number of issues.
#' @param verbose (logical) If `TRUE` (default), displays messages during the checking process.
#' @return (int) Number of issues identified if `return_n_issues` is `TRUE`.
#' @examples {
#' \dontrun{
#' # For Olink data
#' check_crossfile_validation(r_o = results_olink,
#'                            m_s = metadata_samples_olink,
#'                            m_p = metadata_proteins_olink,
#'                            assay_type = "olink")
#'
#' # For LAB data
#' check_crossfile_validation(r_o = results_lab,
#'                            m_s = metadata_samples_lab,
#'                            m_p = metadata_analyte_lab,
#'                            assay_type = "lab")
#' }
#' }
#' @export
check_crossfile_validation <- function(r_o,
                                       m_s,
                                       m_p,
                                       assay_type = c("olink", "lab"),
                                       return_n_issues = FALSE,
                                       verbose = TRUE) {
  
  # Match the assay_type argument
  assay_type <- match.arg(assay_type)
  
  # Initialize issue count
  ic <- 0
  
  # Determine identifier columns based on assay type
  if (assay_type == "olink") {
    analyte_id_col <- "olink_id"
    sample_id_col <- "sample_id"
  } else if (assay_type == "lab") {
    analyte_id_col <- "analyte_name"
    sample_id_col <- "sample_id"
  }
  
  # Validate that required columns are present
  if (!analyte_id_col %in% colnames(r_o)) {
    if (verbose) message("   - (-) `", analyte_id_col, "` column missing in results data frame: FAIL")
    ic <- ic + 1
  }
  if (!sample_id_col %in% colnames(m_s)) {
    if (verbose) message("   - (-) `", sample_id_col, "` column missing in metadata sample data frame: FAIL")
    ic <- ic + 1
  }
  if (!analyte_id_col %in% colnames(m_p)) {
    if (verbose) message("   - (-) `", analyte_id_col, "` column missing in metadata analyte data frame: FAIL")
    ic <- ic + 1
  }
  
  # Proceed only if essential columns are present
  if (ic == 0) {
    # Validate sample IDs
    results_sample_ids <- setdiff(names(r_o), analyte_id_col)
    metadata_sample_ids <- unique(m_s[[sample_id_col]])
    
    # Check for discrepancies in sample IDs
    samples_in_results_not_in_metadata <- setdiff(results_sample_ids, metadata_sample_ids)
    samples_in_metadata_not_in_results <- setdiff(metadata_sample_ids, results_sample_ids)
    
    if (length(samples_in_results_not_in_metadata) > 0) {
      if (verbose) message("   - (-) Samples in results missing from metadata samples: FAIL")
      if (verbose) message("\t\t - ", paste(samples_in_results_not_in_metadata, collapse = "\n\t\t - "))
      ic <- ic + 1
    }
    if (length(samples_in_metadata_not_in_results) > 0) {
      if (verbose) message("   - (-) Samples in metadata samples missing from results: FAIL")
      if (verbose) message("\t\t - ", paste(samples_in_metadata_not_in_results, collapse = "\n\t\t - "))
      ic <- ic + 1
    }
    if (length(samples_in_results_not_in_metadata) == 0 && length(samples_in_metadata_not_in_results) == 0) {
      if (verbose) message("  + (+) All sample IDs match between results and metadata samples: OK")
    }
    
    # Validate analyte IDs
    results_analyte_ids <- r_o[[analyte_id_col]]
    metadata_analyte_ids <- m_p[[analyte_id_col]]
    
    analytes_in_results_not_in_metadata <- setdiff(results_analyte_ids, metadata_analyte_ids)
    analytes_in_metadata_not_in_results <- setdiff(metadata_analyte_ids, results_analyte_ids)
    
    if (length(analytes_in_results_not_in_metadata) > 0) {
      if (verbose) message("   - (-) Analytes in results missing from metadata analytes: FAIL")
      if (verbose) message("\t\t - ", paste(analytes_in_results_not_in_metadata, collapse = "\n\t\t - "))
      ic <- ic + 1
    }
    if (length(analytes_in_metadata_not_in_results) > 0) {
      if (verbose) message("   - (-) Analytes in metadata analytes missing from results: FAIL")
      if (verbose) message("\t\t - ", paste(analytes_in_metadata_not_in_results, collapse = "\n\t\t - "))
      ic <- ic + 1
    }
    if (length(analytes_in_results_not_in_metadata) == 0 && length(analytes_in_metadata_not_in_results) == 0) {
      if (verbose) message("  + (+) All analyte IDs match between results and metadata analytes: OK")
    }
  }
  
  if (return_n_issues) return(ic)
} # End of check_crossfile_validation



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title check failed samples file for reported missing vial label ids
#'
#' @description check failed samples file for reported missing vial label ids
#' @param input_results_folder (char) input path folder
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (vector) failed reported ids
#' @export
check_failedsamples <- function(input_results_folder,
                                verbose = TRUE){
  
  filepattern <- "metadata_failedsamples.*.txt"
  
  # Get file matching pattern
  file_metametabolites <- list.files(input_results_folder,
                                     pattern=filepattern,
                                     full.names=TRUE,
                                     recursive = TRUE,
                                     ignore.case = TRUE)
  
  # Check if file is found and deal with many files
  if(length(file_metametabolites) != 1){
    if(length(file_metametabolites) >= 1){
      if(verbose) message("   - (-) `open_file`: more than one file detected: FAIL")
      if(verbose) message("\n     - ", paste0("`",file_metametabolites,"`", collapse = "\n     - "))
    }else{
      if(verbose) message("   + ( ) File [`", filepattern, "`] not found")
      if(verbose) message("   + ( ) NO FAILED SAMPLES reported")
    }
    flag <- FALSE
    return(NULL)
  }else{
    flag <- TRUE
    ofile <- read.delim(file_metametabolites[1], stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  if(flag){
    if(nrow(ofile) == 0){
      if(verbose) message("   + ( ) NO FAILED SAMPLES reported")
      return(NULL)
    }else{
      if("sample_id" %in% colnames(ofile)){
        if(verbose) message("   + ( ) Failed samples reported:\n\t - ", paste(ofile$sample_id, collapse = "\n\t - ") )
        return(as.character(ofile$sample_id))
      }else{
        if(verbose) message("   - (-) `sample_id` column not found: FAIL")
        return(NULL)
      }
    }
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title check metadata phase file
#'
#' @description check the existence of the metadata phase file
#' @param input_results_folder (char) input path folder
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (logical) `TRUE` exist, `FALSE` does not
#' @export
#' 
check_metadata_phase_file <- function(input_results_folder, 
                                      verbose){
  
  batch <- get_full_path2batch(input_results_folder)
  
  file_phase <- list.files(normalizePath(batch),
                           pattern="metadata_phase.txt",
                           ignore.case = TRUE,
                           full.names=TRUE,
                           recursive = TRUE)
  
  # To be adjusted if two different batches are provided:
  if ( (purrr::is_empty(file_phase)) ){
    if(verbose) message("- (-) `BATCH#_YYYYMMDD/metadata_phase.txt` file does not exist: FAIL")
    return(FALSE)
  }else{
    if(length(file_phase) > 1){
      if(verbose) message("- (-) `More than one `metadata_phase.txt` file available. Only one is valid (place the valid one in the BATCH folder): FAIL")
      return(FALSE)
    }else{
      return(TRUE) 
    }
  }
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Check for Missing Values in a Data Frame Column
#'
#' This function checks for missing values in a specified column of a data frame.
#' It returns TRUE if there are no missing values in the column, and FALSE otherwise.
#'
#' @param df A data frame in which the column to be checked is located.
#' @param column The name of the column to check for missing values, as a string.
#'
#' @return A boolean value; TRUE if the specified column has no missing values, FALSE if it does.
#'
#' @examples
#' data <- data.frame(a = c(1, 2, NA, 4), b = c("A", "B", "C", "D"))
#' check_missing_values(data, "a") # returns TRUE
#' check_missing_values(data, "b") # returns FALSE
#'
#' @export
check_missing_values <- function(df, column) {
  if (!is.data.frame(df)) {
    stop("The first argument must be a data frame.")
  }
  
  if (!column %in% names(df)) {
    stop("The specified column does not exist in the data frame.")
  }
  
  # Check for missing values
  return(any(is.na(df[[column]])))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Validate vial labels from DMAQC
#'
#' @description validate vial label from DMAQC
#' @param vl_submitted (vector) results df
#' @param dmaqc_shipping_info (data.frame) dmaqc shipping information
#' @param tissue_code (char) tissue code
#' @param cas (char) CAS code
#' @param phase (char) phase code
#' @param failed_samples (char) metadata_metabolites df
#' @param return_n_issues (logical) if `TRUE` returns the number of issues
#' @param out_qc_folder (char) output qc folder (it creates the folder if it doesn't exist)
#' @param outfile_missed_viallabels (char) file name for missed vial labels
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @export
check_viallabel_dmaqc <- function(vl_submitted,
                                  dmaqc_shipping_info,
                                  tissue_code,
                                  cas,
                                  phase,
                                  failed_samples,
                                  out_qc_folder = NULL,
                                  outfile_missed_viallabels,
                                  return_n_issues = FALSE,
                                  verbose = TRUE){
  
  # issue_count
  ic <- NA
  
  # Remove redundant samples
  vl_submitted <- unique(vl_submitted[!grepl("\\.", vl_submitted)])
  
  # Adjustment for Gerszten lab
  
  if(tolower(cas) == "broad_rg"){
    cas <- "broad"
  }
  
  # There might be multiple phases to check: load both
  ph <- unlist(strsplit(phase, split = "\\|"))
  dmaqc_labels <- vector()
  month <- NULL
  tr <- NULL
  sp <- NULL
  for(i in 1:length(ph)){
    eph <- ph[i]
    if (grepl("HUMAN", eph)) {
      if (lengths(gregexpr("-", eph)) == 2) {
        # this must be HUMAN-MAIN-TR02 FORMAT
        sp <- gsub("(.*)(-)(.*)(-)(.*)", "\\1", eph)
        tr <- gsub("(.*)(-)(.*)(-)(.*)", "\\5", eph)
      } else if (lengths(gregexpr("-", eph)) == 1) {
        # HUMAN PRECOVID
        sp <- gsub("(.*)(-)(.*)", "\\1", eph)
        tr <- "00"
      } else{
        warning(paste("PROBLEM WITH THE PHASE:", eph, ": it doesn't contain a valid format for HUMAN phase"))
      }
    } else if (grepl("PASS", eph)) {
      sp <- gsub("(.*)(-)(.*)", "\\1", eph)
      month <- gsub("(.*)(-)(.*)", "\\3", eph)
      month <- as.integer(month)
    }

    dmaqc_shipping_df <- read.delim(dmaqc_shipping_info, stringsAsFactors = FALSE)
    
    if(tolower(sp) == "human"){
      dmaqc_labels_temp <- dmaqc_shipping_df$vial_label[which(
        dmaqc_shipping_df$bic_tissue_code == tissue_code &
          dmaqc_shipping_df$site_code == tolower(cas) &
          dmaqc_shipping_df$phase == sp & 
          dmaqc_shipping_df$tranche == tr
      )]
    }else{
      dmaqc_labels_temp <- dmaqc_shipping_df$vial_label[which(
        dmaqc_shipping_df$bic_tissue_code == tissue_code &
          dmaqc_shipping_df$site_code == tolower(cas) &
          dmaqc_shipping_df$phase == sp &
          dmaqc_shipping_df$animal_age == month
      )]
    }

    if(i == 1){
      dmaqc_labels <- as.character(dmaqc_labels_temp)
    }else{
      dmaqc_labels <- append(dmaqc_labels, as.character(dmaqc_labels_temp))
    }
  }
  
  if( length(dmaqc_labels) == 0){
    if(verbose) message("  + (+) DMAQC CHECK POINT: sample IDs not available in DMAQC dataset. Most frequent cause of the error: 
                        - Does the tissue code for this folder structure contain the right tissue code?
                        - Are you using the right phase code? E.g., 'human', or 'pass1a-06'
                        - Did you provide the right code for the cas site? (for example, `broad`, instead of `broad_prot`)
                        
                        Otherwise, it needs to be revised with DMAQC")
    ic <- "NOT_AVAILABLE"
  }else{
    # Report summary counts for clarity
    n_expected <- length(dmaqc_labels)
    n_submitted <- length(vl_submitted)
    n_matched <- length(intersect(vl_submitted, dmaqc_labels))
    if(verbose) message("   + ( ) DMAQC summary: ", n_expected, " expected, ", n_submitted, " submitted, ", n_matched, " matched")
    
    if( setequal(vl_submitted, dmaqc_labels) ){
      if(verbose) message("  + (+) DMAQC CHECK POINT: samples sent to CAS have been processed: OK")
      ic <- "OK"
    }else{
      # Initialize ic as OK, will be set to FAIL if issues found
      ic <- "OK"
      
      # CHECK: samples expected by DMAQC but not in submission
      samples_missed <- setdiff(dmaqc_labels, vl_submitted)
      if( !(is.null(failed_samples) & purrr::is_empty(samples_missed)) ) {
        if( all(samples_missed %in% failed_samples) ){
          if(verbose) message("  + (+) DMAQC CHECK POINT: samples sent to CAS have been processed (with known issues for some samples): OK")
        }else{
          samplesmissedonly <- samples_missed[!(samples_missed %in% failed_samples)]
          if(verbose){
            message("   - (-) DMAQC CHECK POINT: samples not found in `metadata_results`: FAIL")
            message("\t - ", paste(samplesmissedonly, collapse = "\n\t - "))
          }
          missed_out <- data.frame(vial_label = samplesmissedonly)
          missed_out$cas <- cas
          
          # Create output folder
          if (is.null(out_qc_folder)){
            out_qc_folder <- getwd()
          }else{
            out_qc_folder <- create_folder(out_qc_folder)
          }
          
          out_plot_large <- file.path(normalizePath(out_qc_folder), paste0(outfile_missed_viallabels,"-missed_viallabels-in-cas.txt"))
          write.table(missed_out, out_plot_large, row.names = FALSE, sep = "\t", quote = FALSE)
          if(verbose) message("   - ( ) File `", paste0(outfile_missed_viallabels,"-missed_viallabels-in-cas.txt"), "` available with missed vial labels")
          ic <- "FAIL"  
        }
      }else{
        if( !purrr::is_empty(samples_missed) ){
          if(verbose){
            message("   - (-) DMAQC CHECK POINT: samples not found in `metadata_results`: FAIL")
            message("\t - ", paste(samples_missed, collapse = "\n\t - "))
          }
          missed_out <- data.frame(vial_label = samples_missed)
          missed_out$cas <- cas
          out_plot_large <- file.path(normalizePath(out_qc_folder), paste0(outfile_missed_viallabels,"-missed_viallabels-in-cas.txt"))
          write.table(missed_out, out_plot_large, row.names = FALSE, sep = "\t", quote = FALSE)
          if(verbose) message("   - ( ) File ", paste0(outfile_missed_viallabels,"-missed_viallabels-in-cas.txt"), " available with missed vial labels")
          ic <- "FAIL"
        }
      }
      # CHECK: extra samples coming in a submission (not available in DMAQC)
      samples_extra <- setdiff(vl_submitted, dmaqc_labels)
      if(!purrr::is_empty(samples_extra)){
        if(verbose){
          message("   - (-) DMAQC CHECK POINT: CAS SITE IS PROVIDING ", length(samples_extra), " SAMPLE IDS THAT ARE NOT IN DMAQC: REVISE!")
          message("\t - ", paste(samples_extra, collapse = "\n\t - "))
        }
        ic <- "FAIL"
      }
    }
  }
  
  if(return_n_issues) return(ic)
}


#' @title extract ASSAY from input folder path
#'
#' @description extract ASSAY from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) ASSAY code
#' @export
validate_assay <- function(input_results_folder){
  
  assay <- stringr::str_extract(string = input_results_folder,
                                pattern = "(?<=T\\d{2}/)(IONPNEG|RPNEG|RPPOS|HILICPOS|LRPPOS|LRPNEG|3HIB|AA|AC_DUKE|ACOA|BAIBA|CER_DUKE|KA|NUC|OA|SPHM|OXYLIPNEG|ETAMIDPOS|AC_MAYO|AMINES|CER_MAYO|TCA|LAB_GLC|LAB_INS|PROT_PH|PROT_PR|PROT_AC|PROT_UB|PROT_OL|PROT_OX|LAB_CK|LAB_CRT|LAB_CONV)")
  if(is.na(assay)){
    stop("ASSAY not found in the folder structure")
  }else{
    return(assay)
  }
}

#' @title extract BATCH_YYYYMMDD folder
#'
#' @description extract BATCH_YYYYMMDD folder from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) BATCH_YYYYMMDD folder name
#' @export
validate_batch <- function(input_results_folder){
  
  batch_folder <- stringr::str_extract(string = input_results_folder, 
                                       pattern = "(.*/BATCH\\d{1,2}\\_\\d{8})/")
  
  if(is.na(batch_folder)){
    stop("`BATCH#_YYYYMMDD` folder is not recognized in the folder structure.")
  }else{
    return(batch_folder)
  }
}


#' @title validate cas code
#'
#' @description validate CAS code
#' @param cas (char) cas code
#' @export
validate_cas <- function(cas){
  valid_cas_sites <- c("mssm",
                       "broad_met",
                       "bic",
                       "emory",
                       "mayo",
                       "stanford",
                       "umichigan",
                       "gtech",
                       "duke",
                       "pnnl",
                       "broad_prot",
                       "broad_rg")
  if(!(cas %in% valid_cas_sites)){
    stop("cas: <", cas, "> is not valid. Must be one of the following:\n - ", paste(valid_cas_sites, collapse = "\n - "))
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Validate date-time format in a data frame column
#' 
#' @description
#' This function validates that all the values in a specified column of a given data frame
#' adhere to the date-time format: "MM/DD/YYYY HH:MM:SS AM/PM". If any value does not comply 
#' with this format, it prints out those values.
#' 
#' @param df A data frame containing the column to be validated.
#' @param column_name The name of the column in `df` which contains the date-time values.
#' @param verbose Logical. If TRUE, messages are printed to the console.
#' 
#' @return 
#' This function returns the number of issues detected
#' 
#' @examples 
#' 
#' df <- data.frame(id = 1:6,
#'                  datetime = c("12/31/2023 11:59:59 PM", "01/01/2024 00:00:00 AM", 
#'                                "02/29/2024 12:00:00 PM", "13/01/2024 01:00:00 PM", 
#'                                "02/28/2024 24:00:00 PM", "02/30/2024 12:00:00 PM"))
#' validate_dates_times(df, "datetime", TRUE)
#' 
#' @export
validate_dates_times <- function(df, column_name, verbose = TRUE) {
  
  # issue count
  ic <- 0
  
  # Check if the column exists in the dataframe
  if (!column_name %in% names(df)) {
    stop(paste("Column", column_name, "does not exist in the data frame."))
  }
  
  # Check NA and empty values:
  icna <- validate_na_empty(df = df, col_name = column_name, verbose = verbose)
  ic <- ic + icna
  
  # Validate the dates with the specified format
  datetime_values <- df[[column_name]]
  parsed_datetimes <- lubridate::parse_date_time(datetime_values, 
                                                 orders = c("mdy HM", 
                                                            "mdy HMp", 
                                                            "mdy H:M:S",
                                                            "m/d/y h:M:s a", 
                                                            "m/d/y h:M a",
                                                            "mdY h:M:S p"), 
                                                 quiet = TRUE)
  
  # Detect incorrect format
  incorrect_format <- is.na(parsed_datetimes)
  
  # Print the incorrect dates
  incorrect_values <- datetime_values[incorrect_format]
  if (length(incorrect_values) > 0) {
    if(verbose) message("   - (-)`", column_name, "`: Values in incorrect format: `", paste(incorrect_values, collapse = ", "), "`")
      ic <- ic + 1
  }else{
    if(verbose) message("  + (+) All dates are valid.")
  }
  
  # Return the result
  return(ic)
}


#' @title Validate 'lc_column_id' column
#'
#' @description
#'  This function checks the 'lc_column_id' column of a provided data frame
#' to ensure that it exists, contains no NA values, and contains no spaces
#' in its entries. It also reports the number of unique values in the column.
#'
#' @param df A data frame that should contain the 'lc_column_id' column.
#' @param column_name The name of the column in `df` which contains the date-time values.
#' @param verbose A logical indicating whether to print informative messages.
#' Default is TRUE.
#'
#' @return An invisible NULL. The function is used mainly for its side effects
#' (i.e., printing validation results).
#'
#' @examples
#' df <- data.frame(lc_column_id = c("id1", "id2", "id3", "id1", "id 2", NA))
#' validate_lc_column_id(df, column_name = "lc_column_id")
#' 
#' @export
validate_lc_column_id <- function(df, 
                                  column_name, 
                                  verbose = TRUE) {
  
  # issue counter
  ic <- 0
  
  # Check if the column exists in the dataframe
  if (!column_name %in% names(df)) {
    stop(paste("Column", column_name, "does not exist in the data frame."))
  }
  
  # Check NA and empty values:
  icna <- validate_na_empty(df = df, col_name = column_name, verbose = verbose)
  ic <- ic + icna
  
  # check for spaces in values
  if (any(grepl(" ", df[[column_name]]))) {
    if(verbose) message("   - (-) Spaces detected in column `", column_name, "`: FAIL")
    ic <- ic + 1
  }
  
  # report number of unique values
  num_unique <- length(unique(df[[column_name]]))
  if(verbose) message(paste("   - ( ) Number of unique values in column `", column_name, "`: ", num_unique))
  if(num_unique > 3){
    if(verbose) message(paste("   - (!) Warning: the number of LC columns might be too high. Please, revise "))
  }
  
  return(ic)
}


#' Validate Column for NA and Empty Values
#'
#' This function checks if a specified column in a data frame contains either NA or empty values.
#'
#' @param df A data frame.
#' @param col_name A character string specifying the name of the column to check.
#' @param verbose A logical indicating whether to print informative messages. Default is TRUE.
#'
#' @return Number of issues
#'
#' @examples
#' df <- data.frame(A = c("a", "", NA, "d"), B = 1:4)
#' validate_na_empty(df, "A")
#'
#' @export
validate_na_empty <- function(df, col_name, verbose = TRUE) {
  
  # issue count
  ic <- 0
  
  # check if col_name is in column names
  if (!col_name %in% colnames(df)) {
    if(verbose) message(paste("   - (-) Column `", col_name, "` not found in the data frame: FAIL"))
    stop("This column ", col_name, " does not exist")
  }
  
  # check for NA values
  if (any(is.na(df[[col_name]]))) {
    if(verbose) message(paste("   - (-) NA values detected in column `", col_name, "`: FAIL"))
    ic <- ic + 1
  }
  
  # check for empty values
  # check for empty values, ignoring NA
  if (any(df[[col_name]][!is.na(df[[col_name]])] == "")) {
    if(verbose) message(paste("   - (-) Empty values detected in column `", col_name, "`: FAIL"))
    ic <- ic + 1
  }
  
  return(ic)
}


#' @title Extract PHASE from input folder path
#'
#' @description extract ASSAY from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @param return_phase (char) return the phase only if `TRUE` (default)
#' @return (vector) PHASE code
#' @export
validate_phase <- function(input_results_folder, return_phase = TRUE){
  phase <- stringr::str_extract(string = input_results_folder,
                                pattern = "(PASS1A-06|PASS1A-18|PASS1B-06|PASS1B-18|PASS1C-06|PASS1C-18|PASS1AC-06|HUMAN|HUMAN-PRECOVID|HUMAN-MAIN)")
  if( is.na(phase) | phase == "NA" ){
    stop("- (-) Project phase is not found in the folder structure. Please, check the MoTrPAC control vocabulary guidelines")
  }else{
    if(return_phase) return(phase)
  }
}


#' @title extract PROCESSED_YYYYMMDD folder from input folder path
#'
#' @description extract PROCESSED_YYYYMMDD folder from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) PROCESSED_YYYYMMDD folder name
#' @export
validate_processFolder <- function(input_results_folder){

  processfolder <- stringr::str_extract(string = input_results_folder, pattern = "PROCESSED_[0-9]+")

  if(is.na(processfolder)){
    processfolder <- stringr::str_extract(string = input_results_folder, pattern = "(BIC){0,1}RESULTS_[0-9]+")
  }

  if(is.na(processfolder)){
    stop("`PROCESS_YYYYMMDD` or `RESULTS_YYYYMMDD` folder is not recognize in the folder structure")
  }else{
    return(processfolder)
  }
}


#' @title extract and validate TISSUE CODE from input folder path
#'
#' @description extract and validate TISSUE CODE from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) PHASE code
#' @export
validate_tissue <- function(input_results_folder){
  tissue_code <- gsub("(.*)(T[0-9]{2,3})(.*)", "\\2", input_results_folder)

  if(!tissue_code %in% bic_animal_tissue_code$bic_tissue_code){
    stop("tissue_code: `", tissue_code, "` is not valid. Must be one of the following codes (check data object `MotrpacBicQC::bic_animal_tissue_code`):\n- ", paste(bic_animal_tissue_code$bic_tissue_code, collapse = ", "))
  }else{
    return(tissue_code)
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title validate a phase with two phases (pass1a and 1c)
#' 
#' @description This function only works to validate two phases reported
#' as for example, 'PASS1A-06|PASS1C-06' using the separator '|'
#' @param phase_details (char) expected output of `set_phase`
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (char) the expected phase_details function
#' @export
validate_two_phases <- function(phase_details,
                                verbose = TRUE){
  
  if( !grepl("\\|", phase_details) ) 
    stop("This function only validate two phases submitted (e.g PASS1A-06|PASS1C-06), i.e., the variable does not contain the separator '|' required to report two phases")
  
  phase1 <- gsub("(PASS1[A|C]\\-\\d{2})(|.*)", "\\1", phase_details)
  validate_phase(phase1, return_phase = FALSE)
  phase2 <- gsub("(.*|)(PASS1[A|C]\\-\\d{2})", "\\2", phase_details)
  validate_phase(phase2, return_phase = FALSE)
  # only pass1a and 1c expected for this case]
  pass1st <- gsub("(PASS1[A|C])(\\-\\d{2}|.*)", "\\1", phase_details)
  age1st <- gsub("(PASS1[A|C]\\-)(\\d{2})(|.*)", "\\2", phase_details)
  pass2nd <- gsub("(.*|)(PASS1[A|C])(\\-\\d{2})", "\\2", phase_details)
  age2nd <- gsub("(.*|)(PASS1[A|C]\\-)(\\d{2})", "\\3", phase_details)
  
  if(age1st != age2nd){
    stop(paste(phase_details), ": the phase ages reported in `metadata_phase.txt` don't match for these 2 phases: MUST BE CORRECTED")
  }
  if(pass1st == pass2nd){
    stop(paste(phase_details), ": the two reported phases in `metadata_phase.txt` are the same: MUST BE CORRECTED")
  }
  if(verbose) return("Two phases reported and they are ok")
}

#' Validate UniProt IDs
#'
#' This function checks if a given vector of IDs are valid UniProt IDs.
#' It removes NA values, empty strings, and white spaces before validation.
#' Each remaining ID is then checked against the UniProt database.
#' The function outputs a boolean value indicating whether all IDs are valid.
#' It also prints messages for any ID that fails the validation.
#'
#' @param ids A character vector of potential UniProt IDs.
#'
#' @return A boolean value. TRUE if all non-NA, non-empty, and non-whitespace
#'         IDs in the input vector are valid UniProt IDs; FALSE otherwise.
#'
#' @examples
#' # VALID
#' ids1 <- c("P12345", "Q67890", NA, "", " ")
#' if(validate_uniprot_ids_with_uniprot(ids1)) print("Valid UNIPROT IDs")
#' ids2 <- c("P12345", "Q67890", NA, "", " ", "iamwrong")
#' if(!validate_uniprot_ids_with_uniprot(ids1)) print("Invalid UNIPROT IDs")
#'
#' @note
#' This function requires an internet connection to access the UniProt database.
#' It uses the 'httr' package for HTTP requests.
#' The function can be slow for large datasets due to multiple web requests.
#' Also, be aware of potential rate limits or access restrictions on the UniProt API.
#'
#' @export
validate_uniprot_ids_with_uniprot <- function(ids) {
  base_url <- "https://www.uniprot.org/uniprot/"
  
  # Remove NAs, empty strings, and strings with only white spaces
  ids <- ids[!is.na(ids) & ids != "" & ids != " "]
  
  # Split concatenated IDs and create a unique list of IDs
  split_ids <- unlist(strsplit(ids, "_"))
  unique_ids <- unique(split_ids)
  
  # Function to check a single ID
  check_id <- function(id) {
    response <- httr::GET(paste0(base_url, id, ".txt"))
    is_uniprot <- httr::status_code(response) == 200
    if(!is_uniprot) message(paste0("\t- (-) UNIPROT ENTRY `",  id, "` NOT VALID: FAIL"))
    return(is_uniprot)
  }
  
  # Apply check_id function to all unique IDs
  all_true <- all(sapply(unique_ids, check_id))
  return(all_true)
}



#' Validate Dates in a Specified Column of a Data Frame
#'
#' This function checks for the validity of dates in a specified column of a given data frame. 
#' Valid dates are in the format YYYY-MM-DD, with year values between 2018 and 2026, 
#' month values between 1 and 12, and day values between 1 and 31. 
#' The function prints a list of invalid dates and a success message if all dates are valid.
#'
#' @title Validate YYYY-MM-DD Dates in a Data Frame
#'
#' @param df A data frame that contains the date information to be validated.
#' @param date_column A character string specifying the name of the column in `df` that contains the dates to be validated.
#' @param verbose A logical value indicating whether or not to print messages (default: `TRUE`).
#'
#' @return number of issues found
#'
#' @examples
#' df <- data.frame(
#'   extraction_date = c("2022-01-31", "2023-12-01", "2025-11-30"),
#'   other_column = 1:3
#' )
#' ic <- validate_yyyymmdd_dates(df, "extraction_date")
#' @export
validate_yyyymmdd_dates <- function(df, date_column, verbose = TRUE) {
  
  # set issue count
  ic <- 0
  
  # Check if date_column exists in df
  if(!(date_column %in% colnames(df))){
    stop(paste0("Column ", date_column, " not found in the data frame."))
    ic <- ic + 1
    return(ic)
  }
  
  # Extract the date column
  date_vector <- df[[date_column]]
  
  # Check NA and empty values:
  icna <- validate_na_empty(df = df, col_name = date_column, verbose = verbose)
  ic <- ic + icna
  
  # Check dash intead of -
  check_dash <- grepl("\\/", date_vector)
  if(any(check_dash)){
    if(verbose) message("   - (-)`", date_column, "`: Invalid dates detected using `/` instead of `-`: ", paste(date_vector[check_dash], collapse = ", "))
    ic <- ic + 1
    return(ic)
  }
  
  # Check for invalid date format
  incorrect_format <- !grepl("^\\d{4}-\\d{2}-\\d{2}$", date_vector)
  
  # Check for invalid year, month, or day
  split_dates <- strsplit(date_vector, "-")
  incorrect_components <- sapply(split_dates, function(date_parts) {
    year <- as.integer(date_parts[1])
    month <- as.integer(date_parts[2])
    day <- as.integer(date_parts[3])
    
    year_out_of_range <- year < 2018 | year > 2026
    month_out_of_range <- month < 1 | month > 12
    day_out_of_range <- day < 1 | day > 31
    
    year_out_of_range | month_out_of_range | day_out_of_range
  })
  
  # Combine results
  incorrect_dates <- incorrect_format | incorrect_components
  
  if(any(incorrect_dates)){
    if(verbose) message("   - (-) `", date_column, "`: Invalid dates detected: ", paste(date_vector[incorrect_dates], collapse = ", "))
    ic <- ic + 1
  } else {
    if(verbose) message("  + (+) `", date_column, "`: All dates are valid.")
  }
  
  return(ic)
}

