
# All about the RefMet data dictionary

#' @title Get and Validate the Entire RefMet Database from Metabolomics Workbench
#'
#' @description This function fetches and validates the Metabolomics Data Dictionary 
#' from the Metabolomics Workbench. It provides options to remove duplicates.
#'
#' @param remove_duplications Logical; if `TRUE`, removes duplicate entries based on 
#' the `refmet_name` column.
#' @param verbose Logical; if `TRUE` (default), displays progress messages and warnings 
#' during the function execution.
#'
#' @return Returns a data frame with the following columns:
#' \describe{
#'   \item{\code{refmet_name}}{Character; the name standarized refmet name}
#'   \item{\code{pubchem_cid}}{Character; the PubChem compound ID.}
#'   \item{\code{lm_id}}{Character; the LIPID MAPS ID.}
#'   \item{\code{inchi_key}}{Character; the International Chemical Identifier Key.}
#'   \item{\code{exactmass}}{Numeric; the exact mass of the metabolite.}
#'   \item{\code{formula}}{Character; the chemical formula of the metabolite.}
#'   \item{\code{super_class}}{Character; the superclass category of the metabolite.}
#'   \item{\code{main_class}}{Character; the main class category of the metabolite.}
#'   \item{\code{sub_class}}{Character; the subclass category of the metabolite.}
#'   \item{\code{hmdb_id}}{Character; the Human Metabolome Database ID.}
#'   \item{\code{kegg_id}}{Character; the Kyoto Encyclopedia of Genes and Genomes ID.}
#' }
#' Each row of the data frame represents a unique metabolite entry from the 
#' Metabolomics Workbench Data Dictionary. 
#'
#' @details This function downloads the entire RefMet database from the Metabolomics 
#' Workbench using their REST API. The data is initially fetched in JSON format and 
#' then converted to a data frame. The function checks for the presence of a 'name' 
#' column in the data frame, renaming it to 'refmet_name' for consistency. It also 
#' provides an option to remove duplicate entries based on the 'refmet_name' column.
#' If duplicates are found and \code{remove_duplications} is `FALSE`, the function will
#' list the duplicated IDs but will not remove them. This can be helpful for reviewing 
#' the data quality and consistency.
#'
#' @examples
#' \dontrun{
#'   refmet <- get_and_validate_mdd(remove_duplications = TRUE, verbose = TRUE)
#'   head(refmet)
#' }
#'
#' @export
get_and_validate_mdd <- function(remove_duplications = FALSE, verbose = TRUE){

  name = NULL
  id = NULL

  if(verbose) message("- Warning: Downloading data from Metabolomics Workbench. This might take a few minutes.")
  # REST metabolomics workbench data dictionary
  # Previous REST versions
  # refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/motrpac")
  # refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/all/")
  refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/all_ids/")
  
  df <- purrr::map_dfr(refmetjson, ~ as.data.frame(.x), .id = "id")

  if( !("name" %in% colnames(df)) ){
    stop("`refmet_name` column not found in the Metabolomics Workbench data dictionary")
  }else{
    df <- rename(df, refmet_name = name)
  }

  # CHECK DUPLICATIONS
  if(any(duplicated(df$refmet_name))){
    duplications <- df[duplicated(df$refmet_name),]
    if(verbose) message("Duplicated ids: ", length(duplications$refmet_name))
    if(verbose) message("IDS: ", paste(duplications$refmet_name, collapse = ", "))
  }

  refmet <- subset(df, select = -c(id))

  if(remove_duplications){
    # REMOVE DUPLICATIONS
    if( any(duplicated(refmet$refmet_name)) ){
      duplirefmets <- length(refmet$refmet_name[(duplicated(refmet$refmet_name))])
      # message("WARNING: [ ",duplirefmets, " ] DUPLICATION(s) found in data dictionary!")
      refmet <- refmet[!(duplicated(refmet$refmet_name)),]
    }
  }
  return(refmet)
}

#' @title Validate refmet_name
#'
#' @description Validate the refmet_name using the Metabolomics Workbench API
#' @param dataf (data.frame) Data frame with the refmet_name column
#' @param verbose (logical) `TRUE` (default) shows messages)
#' @return (numeric) number of refmet_name ids not available in RefMet
#' @export
validate_refmetname <- function(dataf, verbose = TRUE){

  # Input validation: accept a data.frame (or coerce a list with a
  # `refmet_name` element). Without this check, passing a plain list makes
  # `nrow(dataf)` return NULL and the loop silently iterates zero times,
  # returning 0 missed IDs as if everything had validated successfully.
  if(!is.data.frame(dataf)){
    if(is.list(dataf) && "refmet_name" %in% names(dataf)){
      dataf <- as.data.frame(dataf, stringsAsFactors = FALSE)
    }else{
      stop("`dataf` must be a data.frame with a `refmet_name` column ",
           "(got object of class: ", paste(class(dataf), collapse = ", "), ").")
    }
  }
  if(!("refmet_name" %in% colnames(dataf))){
    stop("`dataf` must contain a `refmet_name` column.")
  }
  if(nrow(dataf) == 0){
    if(verbose) message("  + (+) `refmet_name validation`: no rows to validate.")
    return(0)
  }

  irm <- 0
  n_istd <- 0
  n_api_errors <- 0
  n_missing <- 0
  for(i in seq_len(nrow(dataf))){
    rn <- dataf$refmet_name[i]

    # Handle NA / empty names up front. Without this, grepl() below would
    # return NA on missing values and the surrounding `if()` would error
    # with "missing value where TRUE/FALSE needed", aborting the QC run.
    if(is.na(rn) || !nzchar(trimws(rn))){
      if(verbose) message(paste0("   - (-) `refmet_name validation`: row ", i, " has a missing/empty `refmet_name` (Error RN0)"))
      n_missing <- n_missing + 1
      irm <- irm + 1
      next
    }

    # Skip internal standards
    if(grepl("\\[iSTD\\]", rn)){
      n_istd <- n_istd + 1
      next
    }
    
    tosearch <- "_hp_|_rp_|_rn_|_in_|_lp_|_ln_"
    if(grepl(tosearch, rn)){
      rn <- gsub("(.*)(_\\w{2}_\\w{1})", "\\1", rn)
    }
    
    search_api <- paste0("https://www.metabolomicsworkbench.org/rest/refmet/match/",URLencode(rn),"/name/")
    # Don't crash the whole QC run if the Metabolomics Workbench API is
    # unreachable or returns a non-JSON payload (e.g., an HTML error page).
    here <- tryCatch(
      jsonlite::fromJSON(search_api),
      error = function(e) NULL
    )
    if(is.null(here) || is.null(here$refmet_name)){
      if(verbose) message(paste0("   - (-) `refmet_name validation`: [`", rn, "`] could not be validated (API/network error). Skipped (Error RN3)"))
      n_api_errors <- n_api_errors + 1
      irm <- irm + 1
      next
    }
    if(here$refmet_name == "-"){
      if(verbose) message(paste0("   - (-) `refmet_name validation`: [`", rn, "`] not available in RefMet. Please, contact MW/BIC (Error RN1)"))
      irm <- irm + 1
    }else{
      if(here$refmet_name != rn){
        # Names written as "A/B" may be represented by RefMet as one accepted
        # alternative, so accept matches to any slash-delimited component.
        if(grepl("/", rn, fixed = TRUE)){
          rn_alternatives <- trimws(strsplit(rn, "/", fixed = TRUE)[[1]])
          if(here$refmet_name %in% rn_alternatives){
            next
          }
        }
        if(verbose) message(paste0("   - (-) `refmet_name validation`: [`", rn, "`] must be modified to the RefMet Standardized name: \"", here$refmet_name, "\" (Error RN2)"))
        irm <- irm + 1
      }
    }
  }
  if(n_istd > 0){
    if(verbose) message("  + (+) `refmet_name validation`: Internal standards [iSTD] skipped: ", n_istd)
  }
  if(n_missing > 0){
    if(verbose) message("   - (-) `refmet_name validation`: missing/empty values: ", n_missing)
  }
  if(n_api_errors > 0){
    if(verbose) message("   - (-) `refmet_name validation`: API/network errors (unverified): ", n_api_errors)
  }
  if(irm > 0){
    if(verbose) message("   - (-) `refmet_name validation`: Total number of missed ids on MW: ", irm)
  }
  return(irm)
}
