
# All about the RefMet data dictionary


#' @title Get and validate the entire RefMet database from Metabolomics Workbench
#'
#' @description Get and validate Metabolomics Data Dictionary from
#' Metabolomics Workbench
#' @param remove_duplications (logical) if `TRUE``, removes duplications.
#' @return (vector) PHASE code
#' @export

# get_and_validate_mdd <- function(remove_duplications = FALSE){
#   
#   refmet <- MotrpacBicQC::metabolomics_data_dictionary
# 
#   if(remove_duplications){
#     # REMOVE DUPLICATIONS
#     if( any(duplicated(refmet$refmet_name)) ){
#       duplirefmets <- length(refmet$refmet_name[(duplicated(refmet$refmet_name))])
#       # message("WARNING: [ ",duplirefmets, " ] DUPLICATION(s) found in data dictionary!")
#       refmet <- refmet[!(duplicated(refmet$refmet_name)),]
#     }
#   }
#   return(refmet)
# }

get_and_validate_mdd <- function(remove_duplications = FALSE){

  .id = name = NULL

  # REST metabolomics workbench data dictionary
  # Previous REST version (motrpac only)
  # refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/motrpac")
  refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/all/")

  dt_list <- purrr::map(refmetjson, as.data.table)
  dt <- data.table::rbindlist(dt_list, fill = TRUE, idcol = T)
  df <- as.data.frame(dt)

  colnames(df) <- tolower(colnames(df))

  if( !("name" %in% colnames(df)) ){
    stop("`refmet_name` column not found in the Metabolomics Workbench data dictionary")
  }else{
    df <- rename(df, refmet_name = name)
  }

  # CHECK DUPLICATIONS
  if(any(duplicated(df$refmet_name))){
    duplications <- df[duplicated(df$refmet_name),]
    message("Duplicated ids: ", length(duplications$refmet_name))
    message("IDS: ", paste(duplications$refmet_name, collapse = ", "))
    message("DUPLICATIONS IN REFMET ONLINE (REST VERSION)")
  }

  refmet <- subset(df, select = -c(.id))

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
#' @description Validate the refment_name using the Metabolomics Workbench API
#' @param dataf (data.frame) Data frame with the refmet_name column
#' @param verbose (logical) `TRUE` (default) shows messages)
#' @return (numeric) number of refmet_name ids not available in RefMet
#' @export
validate_refmetname <- function(dataf, verbose){

  irm <- 0
  idna <- 0
  for(i in 1:dim(dataf)[1]){
    rn <- dataf$refmet_name[i]
    
    tosearch <- "_hp_|_rp_|_rn_|_in_|_lp_|_ln_"
    if(grepl(tosearch, rn)){
      # message("\tMulti-isotope detected: ", rn, appendLF = FALSE)
      rn <- gsub("(.*)(_\\w{2}_\\w{1})", "\\1", rn)
      # message(" > ", rn)
    }
    # Remove the standard label
    standards <- "\\[.*\\]"
    if(grepl(standards, rn)){
      # message("\tStandard detected: ", rn, appendLF = FALSE)
      rn <- gsub("(.*)(\\s+\\[iSTD\\])", "\\1", rn)
      # message(" > ", rn)
    }
    
    search_api <- paste0("https://www.metabolomicsworkbench.org/rest/refmet/match/",URLencode(rn),"/name/")
    here <- jsonlite::fromJSON(search_api)
    if(length(here) == 0){
      if(verbose) message(paste0("      (-) `refmet_name` [`", rn, "`] not available in RefMet. Please, contact MW/BIC (Error RN1)"))
      irm <- irm + 1
      idna <- idna + 1
    }
  }
  if(idna > 0){
    if(verbose) message("      (-) Total number of missed ids on MW: ", idna)
  }
  return(irm)
}
