#' @import knitr
#' @import utils

# Remove empty columns in data frame

# ------------------------------------------------------------------------------
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
  file_metametabolites <- list.files(input_results_folder,
                                     pattern=filepattern,
                                     ignore.case = TRUE,
                                     full.names=TRUE,
                                     recursive = TRUE)

  # Check if file is found and deal with many files
  if(length(file_metametabolites) != 1){
    if(length(file_metametabolites) >= 1){
      if(verbose) message("      - (-) PROBLEM: more than one file detected: FAIL")
      if(verbose) message("\n\t\t - ", paste(file_metametabolites, collapse = "\n\t\t - "))
    }else{
      if(verbose) message("      - (-) PROBLEM file [", filepattern, "] not found: FAIL")
    }
    flag <- FALSE
    ofile <- NULL
  }else{
    flag <- TRUE
    ofile <- read.delim(file_metametabolites[1], stringsAsFactors = FALSE, check.names = FALSE)
    ofile <- remove_empty_columns(ofile)
    ofile <- remove_empty_rows(ofile)
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

  list_back <- list("flag" = flag, "df" = ofile)
  return(list_back)
}

# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
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
