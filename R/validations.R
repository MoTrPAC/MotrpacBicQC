# VALIDATIONS

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
    return(TRUE)
  }
  
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
  
  # There might be multiple phases to check: load both
  ph <- unlist(strsplit(phase, split = "\\|"))
  dmaqc_labels <- vector()
  for(i in 1:length(ph)){
    eph <- ph[i]
    pass <- gsub("(.*)(-)(.*)", "\\1", eph)
    if(tolower(pass) != "human"){
      month <- gsub("(.*)(-)(.*)", "\\3", eph)
      month <- as.integer(month)
    }
    
    dmaqc_shipping_df <- read.delim(dmaqc_shipping_info, stringsAsFactors = FALSE)
    
    if(tolower(pass) == "human"){
      dmaqc_labels_temp <- dmaqc_shipping_df$vial_label[which(dmaqc_shipping_df$bic_tissue_code == tissue_code &
                                                                dmaqc_shipping_df$site_code == tolower(cas) &
                                                                dmaqc_shipping_df$phase == pass)]
    }else{
      dmaqc_labels_temp <- dmaqc_shipping_df$vial_label[which(dmaqc_shipping_df$bic_tissue_code == tissue_code &
                                                                dmaqc_shipping_df$site_code == tolower(cas) &
                                                                dmaqc_shipping_df$phase == pass &
                                                                dmaqc_shipping_df$animal_age == month)]
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
    if( setequal(vl_submitted, dmaqc_labels) ){
      if(verbose) message("  + (+) DMAQC CHECK POINT: samples sent to CAS have been processed: OK")
      ic <- "OK"
    }else{
      # CHECK 
      samples_missed <- setdiff(dmaqc_labels, vl_submitted)
      if( !(is.null(failed_samples) & purrr::is_empty(samples_missed)) ) {
        if( all(samples_missed %in% failed_samples) ){
          if(verbose) message("  + (+) DMAQC CHECK POINT: samples sent to CAS have been processed (with known issues for some samples): OK")
          ic <- "OK"
        }else{
          samplesmissedonly <- samples_missed[!(samples_missed %in% failed_samples)]
          if(verbose){
            message("   - (-) DMAQC CHECK POINT: samples not found in `metadata_results`: FAIL")
            message("\t - ", paste(samplesmissedonly, collapse = "\n\t - "))
          }
          missed_out <- data.frame(vial_label = samplesmissedonly)
          missed_out$cas <- cas
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
          message("   - (-) DMAQC CHECK POINT: CAS SITE IS PROVIDING SAMPLES IDS THAT ARE NOT IN DMAQC: REVISE!")
          message("\t - ", paste(samples_extra, collapse = "\n\t - "))
        }        
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
                                pattern = "(IONPNEG|RPNEG|RPPOS|HILICPOS|LRPPOS|LRPNEG|3HIB|AA|AC_DUKE|ACOA|BAIBA|CER_DUKE|CONV|KA|NUC|OA|SPHM|OXYLIPNEG|ETAMIDPOS|AC_MAYO|AMINES|CER_MAYO|TCA|IMM_CRT|IMM_GLC|IMM_INS|PROT_PH|PROT_PR|PROT_AC|PROT_UB)")
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
                       "broad_prot")
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
  parsed_datetimes <- lubridate::parse_date_time(datetime_values, orders = c("mdy HM", "mdy HMp"), quiet = TRUE)
  
  # Detect incorrect format
  incorrect_format <- is.na(parsed_datetimes)
  
  if(verbose){
    # Print the incorrect dates
    incorrect_values <- datetime_values[incorrect_format]
    if (length(incorrect_values) > 0) {
      message("   - (-)`", column_name, "`: Values in incorrect format: `", paste(incorrect_values, collapse = ", "), "`")
      ic <- ic + 1
    }else{
      message("  + (+) All dates are valid.")
    }
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
validate_lc_column_id <- function(df, column_name, verbose = TRUE) {
  
  # issue counter
  ic <- 0
  
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
    return(invisible())
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
  if(is.na(phase)){
    stop("Project phase (e.g. PASS1A-06) is not found in the folder structure, please, check guidelines")
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
#'   extraction_date = c("2022-01-31", "2023-12-01", "2027-11-30"),
#'   other_column = 1:3
#' )
#' validate_yyyymmdd_dates(df, "extraction_date")
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
    if(verbose) message("  + (+) All dates are valid.")
  }
  
  return(ic)
}

