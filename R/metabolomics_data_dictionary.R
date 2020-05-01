
# All about the RefMet data dictionary


#' @title Get and validate Metabolomics Data Dictionary from Metabolomics Workbench
#'
#' @description Get and validate Metabolomics Data Dictionary from
#' Metabolomics Workbench
#' @param remove_duplications (logical) if `TRUE``, removes duplications.
#' @return (vector) PHASE code
#' @export
get_and_validate_mdd <- function(remove_duplications = FALSE){

  .id = NULL

  # REST metabolomics workbench data dictionary
  refmetjson <- jsonlite::fromJSON("https://www.metabolomicsworkbench.org/rest/refmet/motrpac/")

  dt_list <- purrr::map(refmetjson, as.data.table)
  dt <- data.table::rbindlist(dt_list, fill = TRUE, idcol = T)
  df <- as.data.frame(dt)

  colnames(df) <- tolower(colnames(df))

  if( !("refmet_name" %in% colnames(df)) ){
    stop("<refmet_name> column not found in the Metabolomics Workbench data dictionary")
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
