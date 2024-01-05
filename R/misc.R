# Add package by alphabetical order

#' @importFrom data.table rbindlist as.data.table fread
#' @import dplyr
#' @import forcats
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom httr status_code GET
#' @importFrom inspectdf inspect_na
#' @importFrom jsonlite fromJSON
#' @import knitr
#' @importFrom lubridate parse_date_time
#' @import naniar
#' @import progress
#' @import purrr
#' @importFrom readr read_lines
#' @importFrom scales percent
#' @importFrom stats median reorder
#' @import stringr
#' @import tidyr
#' @import utils
#' @import viridis
#____________________________________________________________________________

#' @title Create folder
#'
#' @description Create a directory if it doesn't exist. If no argument is provided,
#' it returns the current working directory
#' @param folder_name (chr) folder name
#' @param verbose (logical) `TRUE` shows messages (default `FALSE`)
#' @examples {
#' create_folder(folder_name = NULL)
#' # Or use this one for a real folder:
#' # create_folder(folder_name = "testing")
#' }
#' @export
create_folder <- function(folder_name = NULL,
                          verbose = FALSE){
  if(!is.null(folder_name)){
    if(!dir.exists(file.path(folder_name))){
      dir.create(file.path(folder_name), recursive = TRUE)
      if(verbose) message("+ Folder `", folder_name,"`created")
      return(folder_name)
    }else{
      return(folder_name)
    }
  }else{
    folder_name <- getwd()
    return(folder_name)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title dl_read_gcp: Data Load, Read file from Google Cloud
#'
#' @description
#' Read a single file from Google Cloud Storage (GSC) into a data table
#'
#' @param path (char) GCS path, i.e., starts with "gs://"
#' @param sep (char) column separator to use with [data.table::fread]
#' @param tmpdir (char) scratch directory to download files from GCS
#' @param gsutil_path (char) path to \code{gsutil} on your computer.
#' Can be "gsutil" if \code{gsutil} is in your \code{$PATH}.
#' @param check_first (char) check if file exists in \code{tmpdir} before
#' downloading it. Read in existing file if it exists.
#' Should be set to \code{TRUE} if you are running this function in parallel.
#' @param header (bool) whether input file has a header line
#' @param verbose (logical) `TRUE` shows messages (default `FALSE`)
#' @param ... optional arguments for [data.table::fread]
#'
#' @return a data table
#'
#' @importFrom data.table fread
#'
#' @examples
#' \dontrun{
#' pheno = dl_read_gcp(path = "gs://your-bucket/file.txt")
#' }
#' @export
dl_read_gcp <- function(path,
                        sep = "\t",
                        header = TRUE,
                        tmpdir,
                        gsutil_path = "gsutil",
                        check_first = TRUE,
                        verbose = FALSE,
                        ...){

  if(!dir.exists(tmpdir)){
    dir.create(tmpdir)
    if(verbose) message(paste0("- New folder ", tmpdir, " created successfully"))
  }

  tmpdir <- normalizePath(tmpdir)

  # Check path
  if(!grepl("gs:\\/\\/", path)){
    stop("The path to the bucket is wrong. Valid example: gs://bucket-name/file-name.csv")
  }else{
    new_path <- file.path(tmpdir, basename(path))
  }
  
  # Detect the operating system
  os_name <- Sys.info()["sysname"]
  
  # Default arguments for Mac
  ignore_std_err <- TRUE
  ignore_std_out <- TRUE
  
  # Change default arguments if the OS is Windows
  if (os_name == "Windows") {
    ignore_std_err <- FALSE
    ignore_std_out <- FALSE
  }

  # only download if it doesn't exist to avoid conflicts when running this script in parallel; clear scratch space when you're done
  if(check_first){
    if( !file.exists(new_path) ){
      cmd <- sprintf('%s cp %s %s', gsutil_path, path, tmpdir)
      if(verbose) message(paste0("- Running command ", cmd))
      system(cmd,
             ignore.stdout = ignore_std_out,
             ignore.stderr = ignore_std_err)
      message("- Downloaded file: ", new_path)
    }else{
      if(verbose) message(paste0("- The file <", new_path, "> already exists"))
    }
  }else{
    if(verbose) message(paste("- Downloading file (from GCP) <", basename(path), ">"))
    cmd <- sprintf('%s cp %s %s', gsutil_path, path, tmpdir)
    system(cmd,
           ignore.stdout = ignore_std_out,
           ignore.stderr = ignore_std_err)
    message("- Downloaded file: ", new_path)
  }
  # read in the data as a data.table
  if(file.exists(new_path)){
    dt <- data.table::fread(new_path, sep=sep, header=header,...)
    return(dt)
  }else{
    stop("- Problems loading the file. Possible reason: the file does not exist in the bucket anymore. Please, validate the address. Re-run this command again with `verbose = TRUE`)")
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Generate the phase detail for submissions
#'
#' @description The phase details is as simple as creating a lower case version
#' of the phase. However, in case of PASS1A/1C a new version has to be generated:
#' pass1ac-06
#' This function detects whether there are two phases, and if so,
#' generate the expected version: either pass1ac-06 or pass1ac-18
#' @param phase_metadata (char) expected output of `set_phase`
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (char) the expected phase_details function
#' @export
generate_phase_details <- function(phase_metadata,
                                   verbose = TRUE){

  if( grepl("\\|", phase_metadata) ){
    pass1st <- gsub("(.*)(\\|.*)", "\\1", phase_metadata)
    animalage <- gsub("(PASS1A\\-)(\\d+)", "\\2", pass1st)
    phase_details <- paste0("pass1ac-", animalage)
  }else{
    phase_details <- tolower(phase_metadata)
  }
  return(phase_details)
}


#' @title Get full path to the batch folder
#'
#' @description Get the full path to the batch folder
#' @param input_results_folder (char) path to the PROCESSED/RESULTS folder to check
#' @return (char) Full path to the `BATCH#_YYYYMMDD` folder
#' @export
get_full_path2batch <- function(input_results_folder){

  batch <- NULL

  if( grepl("(BIC){0,1}RESULTS", input_results_folder) ){
    batch <- gsub("(.*/)((BIC){0,1}RESULTS.*)", "\\1", input_results_folder)
  }else if( grepl("PROCESSED", input_results_folder)){
    batch <- gsub("(.*)(PROCESSED.*)", "\\1", input_results_folder)
  }else{
    stop("   - (-) ERROR: the input results folder missed the PROCESSED or RESULTS folder!")
  }

  return(batch)

}


#' @title filter required columns only
#'
#' @description it returns a data frame with only the required columns for metabolomics and proteomics
#' @param df (data.frame) metadata_metabolites
#' @param type (char) Type of file to filter columns:
#' - `m_m`: metadata metabolites
#' - `m_s`: metadata samples
#' - `v_m`: proteomics vial_metadata
#' - `olproteins`: olink metadata proteins
#' - `olsamples`: olink metadata samples
#' @param name_id (char) specify whether `named` or `unnamed` files
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (data.frame) filtered data frame with only the required columns
#' @examples {
#' df_filtered <- filter_required_columns(df = metadata_metabolites_named, name_id = "named")
#' }
#' @export
filter_required_columns <- function(df,
                                    type = c("m_m", 
                                             "m_s", 
                                             "v_m", 
                                             "olproteins", 
                                             "olsamples"),
                                    name_id = NULL,
                                    verbose = TRUE){

  type <- match.arg(type)

  if (type == "m_m"){
    # Define required columns
    if(name_id == "named"){
      emeta_metabo_coln_named <- c("metabolite_name", "refmet_name", "rt", "mz", "neutral_mass", "formula")
    }else if(name_id == "unnamed"){
      if("neutral_mass" %in% colnames(df)){
        emeta_metabo_coln_named <- c("metabolite_name", "rt", "mz", "neutral_mass")
      }else{
        emeta_metabo_coln_named <- c("metabolite_name", "rt", "mz")
      }
    }else{
      stop("{`name_id`} option not valid. Options: named/unnamed")
    }

    # Now check if present
    colnames(df) <- tolower(colnames(df))
    missing_cols <- setdiff(emeta_metabo_coln_named, colnames(df))
    if (length(missing_cols) > 0) {
      if(verbose) message("   - (-) `metadata_metabolite`: Expected COLUMN NAMES are missed: FAIL")
      message(paste0("\t The following required columns are not present: `", paste(missing_cols, collapse = ", "), "`"))
    } else {
      if(verbose) message("  + (+) All required columns present")
      df <- subset(df, select = emeta_metabo_coln_named)
    }
    return(df)

  } else if (type == "m_s"){
    emeta_sample_coln <- c("sample_id", "sample_type", "sample_order", "raw_file", "extraction_date", "acquisition_date", "lc_column_id")
    missing_cols <- setdiff(emeta_sample_coln, colnames(df))

    if (length(missing_cols) > 0) {
      if(verbose) message("   - (-) `metadata_sample`: Expected COLUMN NAMES are missed: FAIL")
      message(paste0("\t The following required columns are not present: `", paste(missing_cols, collapse = ", "), "`"))
    } else {
      if(verbose) message("  + (+) All required columns present")
      df <- subset(df, select = emeta_sample_coln)
    }
    return(df)

  } else if (type == "v_m"){
    emeta_sample_coln <- c("vial_label", "tmt_plex")
    if( all(emeta_sample_coln %in% colnames(df)) ){
      # deal with tmt11 or tmt16
      if("tmt11_channel" %in% colnames(df)){
        emeta_sample_coln <- append(emeta_sample_coln, "tmt11_channel")
        if(verbose) message("  + (+) All required columns present (tmt11 experiment)")
        df <- subset(df, select = emeta_sample_coln)
      }else if("tmt16_channel" %in% colnames(df)){
        emeta_sample_coln <- append(emeta_sample_coln, "tmt16_channel")
        if(verbose) message("  + (+) All required columns present (tmt16 experiment)")
        df <- subset(df, select = emeta_sample_coln)
      }else{
        if(verbose) message("   - (-) Expected COLUMN NAMES are missed: FAIL")
      }
    }else{
      if(verbose) message("   - (-) Expected COLUMN NAMES are missed: FAIL")
    }
    return(df)
  } else if (type == "olproteins"){
    emeta_sample_coln <- c("olink_id", "uniprot_entry", "assay", "missing_freq", "panel_name", "panel_lot_nr", "normalization")
    missing_cols <- setdiff(emeta_sample_coln, colnames(df))
    
    if (length(missing_cols) > 0) {
      if(verbose) message("   - (-) `metadata_proteins`: Expected COLUMN NAMES are missed: FAIL")
      message(paste0("\t The following required columns are not present: `", paste(missing_cols, collapse = ", "), "`"))
    } else {
      if(verbose) message("  + (+) All required columns present")
      df <- subset(df, select = emeta_sample_coln)
    }
    return(df)
  }else if (type == "olsamples"){
    emeta_sample_coln <- c("sample_id", "sample_type", "sample_order", "plate_id")
    missing_cols <- setdiff(emeta_sample_coln, colnames(df))
    
    if (length(missing_cols) > 0) {
      if(verbose) message("   - (-) `metadata_samples`: Expected COLUMN NAMES are missed: FAIL")
      message(paste0("\t The following required columns are not present: `", paste(missing_cols, collapse = ", "), "`"))
    } else {
      if(verbose) message("  + (+) All required columns present")
      df <- subset(df, select = emeta_sample_coln)
    }
    return(df)
  }
}

#' @title open files
#'
#' @description open files and check that they are right
#' @param input_results_folder (char) input path folder
#' @param filepattern (char) regular expression to find a file in the file system
#' provided
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (list) list with data frame and flag about the status
#' @export
open_file <- function(input_results_folder,
                      filepattern,
                      verbose = TRUE){

  if( !dir.exists(input_results_folder) ){
    flag <- FALSE
    ofile <- NULL
    filename <- NULL
    if(verbose) message("   - (-) The folder doesn't exist: FAIL")
    list_back <- list("flag" = flag, "df" = ofile, "filename" = filename)
    return(list_back)
  }

  # Get file matching pattern
  file_metametabolites <- list.files(normalizePath(input_results_folder),
                                     pattern=filepattern,
                                     ignore.case = TRUE,
                                     full.names=TRUE,
                                     recursive = TRUE)

  # Check if file is found and deal with many files
  if(length(file_metametabolites) != 1){
    if(length(file_metametabolites) >= 1){
      if(verbose) message("   - (-) More than one file detected: FAIL")
      if(verbose) message("\t\t - ", paste(file_metametabolites, collapse = "\n\t\t - "))
    }else{
      if(verbose) message("   - (-) File [`", filepattern, "`] not found: FAIL")
    }
    flag <- FALSE
    ofile <- NULL
    filename <- NULL
  }else{
    flag <- TRUE
    filename <- file_metametabolites[1]
    ofile <- read.delim(filename, stringsAsFactors = FALSE, check.names = FALSE)
    ofile <- remove_empty_columns(ofile, verbose = verbose)
    ofile <- remove_empty_rows(ofile, verbose = verbose)
    if(verbose) message("  + (+) File successfully opened")
  }

  if(flag){
    if(nrow(ofile) == 0){
      if(verbose) message("   - (-) File is empty: FAIL")
      flag <- FALSE
      ofile <- NULL
    }else{
      flag <- TRUE
    }
  }

  list_back <- list("flag" = flag, "df" = ofile, "filename" = filename)
  return(list_back)
}


#' @title remove empty columns
#'
#' @description remove empty columns
#' @param df (char) data frame
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (df) df without empty columns
#' @export
remove_empty_columns <- function(df,
                                 verbose = TRUE){
  df[df == ""] <- NA
  before <- dim(df)[2]
  emptycols <- sapply(df, function (x) all(is.na(x)))
  df <- df[!emptycols]
  after <- dim(df)[2]
  if(before != after){
    n_removed <- before - after
    if(verbose) message("   - (-) ", n_removed, " empty columns found and removed")
    if(verbose) message("\t\t+ Before: ", before, " ->  After: ", after)
  }
  return(df)
}


#' @title remove empty rows in data frame
#'
#' @description remove empty rows in data frame
#' @param df (char) data frame
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (df) df without empty columns
#' @export
remove_empty_rows <- function(df,
                              verbose = TRUE){
  # Remove all rows with NAs or white spaces
  # 1. Check empty spaces and make them NAs
  before <- dim(df)[1]
  df[df == ""] <- NA
  # 2. Remove rows that are all NAs
  df <- df[apply(df, 1, function(x) !all(is.na(x))),]
  after <- dim(df)[1]
  if(before != after){
    n_removed <- before - after
    if(verbose) message("   - (-) ", n_removed, " empty ROWS found and remove")
    if(verbose) message("\t\t+ Before: ", before, " ->  After: ", after)
  }
  return(df)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Set the phase to be validated.
#'
#' @description A group might choose to combine two different phases, due to
#' the complications associated with PASS1A/1C. If they choose to combine
#' two phases, the CAS must provide a new file `metadata_phase.txt` with a single
#' line, as for example: `PASS1A-06|PASS1C-06`. This function checks if the
#' file is available, and set that phase as the phases to validate. In summary,
#' the order of preference is:
#' 1. function's argument: dmaqc_phase2validate (if provided in the validation functions)
#' 2. `metadata_phase.txt` file if available in the batch folder.
#' 3. Phase in folder structure
#' @param input_results_folder (char) path to the PROCESSED/RESULTS folder to check
#' @param dmaqc_phase2validate (data.frame) dmaqc shipping information
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) the phase to be validated.
#' @export
set_phase <- function(input_results_folder,
                      dmaqc_phase2validate,
                      verbose = TRUE){

  phase <- validate_phase(input_results_folder)

  # Check metadata_phase.txt file
  batch <- get_full_path2batch(input_results_folder)

  file_phase <- list.files(normalizePath(batch),
                           pattern="metadata_phase.txt",
                           ignore.case = TRUE,
                           full.names=TRUE,
                           recursive = TRUE)

  # To be adjusted if two different batches are provided:
  if ( !(purrr::is_empty(file_phase)) ){
    phase_details <- readr::read_lines(file_phase, n_max = 1)
    if ( !(is.na(phase_details) || phase_details == '') ){
      if(verbose) message("+ Motrpac phase reported: ", phase_details, " (info from metadata_phase.txt available): OK")

      if( grepl("\\|", phase_details) ){
        validate_two_phases(phase_details = phase_details, verbose = FALSE)
      }

      # And once is checked, proceed...
      if( isFALSE(dmaqc_phase2validate) ){
        dmaqc_phase2validate <- phase_details
      }
    }else{
      if(verbose) message("+ Motrpac phase: ", phase, " (metadata_phase.txt available but EMPTY): FAIL")
      if( isFALSE(dmaqc_phase2validate) ){
        dmaqc_phase2validate <- phase
      }
    }
  }else{
    if(verbose) message("+ Motrpac phase: ", phase, " (metadata_phase.txt file NOT available): FAIL")
    if( isFALSE(dmaqc_phase2validate) ){
      dmaqc_phase2validate <- phase
    }
  }

  return(dmaqc_phase2validate)
}

