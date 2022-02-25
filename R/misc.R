# Add package by alphabetical order

#' @importFrom data.table rbindlist as.data.table
#' @import dplyr
#' @import forcats
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom inspectdf inspect_na
#' @importFrom jsonlite fromJSON
#' @import knitr
#' @import naniar
#' @import progress
#' @import purrr
#' @importFrom scales percent
#' @importFrom stats median reorder
#' @import stringr
#' @import tidyr
#' @import utils
#____________________________________________________________________________

#' @title filter required columns only
#'
#' @description it returns a data frame with only the required columns for metabolomics and proteomics
#' @param df (data.frame) metadata_metabolites
#' @param type (char) Type of file to filter columns:
#' - `m_m`: metadata metabolites
#' - `m_s`: metadata samples
#' - `v_m`: proteomics vial_metadata
#' @param name_id (char) specify whether `named` or `unnamed` files
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (data.frame) filtered data frame with only the required columns
#' @examples {
#' df_filtered <- filter_required_columns(df = metadata_metabolites_named, name_id = "named")
#' }
#' @export
filter_required_columns <- function(df,
                                    type = c("m_m", "m_s", "v_m"),
                                    name_id = NULL,
                                    verbose = TRUE){

  type <- match.arg(type)

  if (type == "m_m"){
    if(name_id == "named"){
      emeta_metabo_coln_named <- c("metabolite_name", "refmet_name", "rt", "mz", "neutral_mass", "formula")
    }else if(name_id == "unnamed"){
      if("neutral_mass" %in% colnames(df)){
        emeta_metabo_coln_named <- c("metabolite_name", "rt", "mz", "neutral_mass")
      }else{
        emeta_metabo_coln_named <- c("metabolite_name", "rt", "mz")
      }
    }else{
      stop("{name_id} option not valid. Options: named/unnamed")
    }

    colnames(df) <- tolower(colnames(df))

    if(all(emeta_metabo_coln_named %in% colnames(df))){
      if(verbose) message("   + (+) All required columns present")
      df <- subset(df, select = emeta_metabo_coln_named)
    }else{
      if(verbose) message("      - (-) Expected COLUMN NAMES are missed: FAIL")
    }
    return(df)
  } else if (type == "m_s"){
    emeta_sample_coln <- c("sample_id", "sample_type", "sample_order", "raw_file")
    if( all(emeta_sample_coln %in% colnames(df)) ){
      if(verbose) message("   + (+) All required columns present")
      df <- subset(df, select = emeta_sample_coln)
    }else{
      if(verbose) message("      - (-) Expected COLUMN NAMES are missed: FAIL")
    }
    return(df)
  } else if (type == "v_m"){
    emeta_sample_coln <- c("vial_label", "tmt_plex")
    if( all(emeta_sample_coln %in% colnames(df)) ){
      # deal with tmt11 or tmt16
      if("tmt11_channel" %in% colnames(df)){
        emeta_sample_coln <- append(emeta_sample_coln, "tmt11_channel")
        if(verbose) message("   + (+) All required columns present (tmt11 experiment)")
        df <- subset(df, select = emeta_sample_coln)
      }else if("tmt16_channel" %in% colnames(df)){
        emeta_sample_coln <- append(emeta_sample_coln, "tmt16_channel")
        if(verbose) message("   + (+) All required columns present (tmt16 experiment)")
        df <- subset(df, select = emeta_sample_coln)
      }else{
        if(verbose) message("      - (-) Expected COLUMN NAMES are missed: FAIL")
      }
    }else{
      if(verbose) message("      - (-) Expected COLUMN NAMES are missed: FAIL")
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

  # Get file matching pattern
  file_metametabolites <- list.files(normalizePath(input_results_folder),
                                     pattern=filepattern,
                                     ignore.case = TRUE,
                                     full.names=TRUE,
                                     recursive = TRUE)

  # Check if file is found and deal with many files
  if(length(file_metametabolites) != 1){
    if(length(file_metametabolites) >= 1){
      if(verbose) message("      - (-) PROBLEM: more than one file detected: FAIL")
      if(verbose) message("\t\t - ", paste(file_metametabolites, collapse = "\n\t\t - "))
    }else{
      if(verbose) message("      - (-) PROBLEM file [", filepattern, "] not found: FAIL")
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
    if(verbose) message("   + (+) File successfully opened")
  }

  if(flag){
    if(nrow(ofile) == 0){
      if(verbose) message("      - (-) File is empty: FAIL")
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
    if(verbose) message("      - (-) ", n_removed, " empty columns found and removed")
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
    if(verbose) message("      - (-) ", n_removed, " empty ROWS found and remove")
    if(verbose) message("\t\t+ Before: ", before, " ->  After: ", after)
  }
  return(df)
}
