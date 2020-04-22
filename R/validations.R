# VALIDATIONS


#' @title validate cas code
#'
#' @description validate CAS code
#' @param cas (char) cas code
#' @export
validate_cas <- function(cas){
  valid_cas_sites <- c("mssm",
                       "broad_met",
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


#' @title extract PROCESSED_YYYYMMDD folder from input folder path
#'
#' @description extract PROCESSED_YYYYMMDD folder from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) PROCESSED_YYYYMMDD folder name
#' @export
validate_processFolder <- function(input_results_folder){

  processfolder <- stringr::str_extract(string = input_results_folder, pattern = "PROCESSED_[0-9]+")

  if(is.na(processfolder)){
    stop("PROCESS_YYYYMMDD folder is not recognize in the folder structure")
  }else{
    return(processfolder)
  }
}



#' @title extract ASSAY from input folder path
#'
#' @description extract ASSAY from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) ASSAY code
#' @export
validate_assay <- function(input_results_folder){

  assay <- stringr::str_extract(string = input_results_folder,
                                pattern = "(IONPNEG|RPNEG|RPPOS|HILICPOS|LRPPOS|LRPNEG|OXYLIPNEG)")
  if(is.na(assay)){
    stop("ASSAY not found in the folder structure")
  }else{
    return(assay)
  }
}



#' @title extract PHASE from input folder path
#'
#' @description extract ASSAY from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) PHASE code
#' @export
validate_phase <- function(input_results_folder){
  phase <- stringr::str_extract(string = input_results_folder,
                                pattern = "(PASS1A-06|PASS1A-18|PASS1B-06|PASS1B-18)")
  if(is.na(phase)){
    stop("<project phase (e.g. PASS1A-06> is not found in the folder structure")
  }else{
    return(phase)
  }
}



#' @title extract and validate TISSUE CODE from input folder path
#'
#' @description extract and validate TISSUE CODE from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @return (vector) PHASE code
#' @export
validate_tissue <- function(input_results_folder){
  tissue_db <- c("T30", "T31", "T32", "T33", "T34", "T35", "T36", "T37",
                 "T38", "T39", "T40", "T41", "T42", "T43", "T44", "T45",
                 "T46", "T47", "T48", "T49", "T50", "T51", "T52", "T53",
                 "T54", "T55", "T56", "T57", "T58", "T59", "T60", "T61",
                 "T62", "T63", "T64", "T65", "T66", "T67", "T68", "T69",
                 "T70")

  tissue_code <- gsub("(.*)(T[0-9]{2})(.*)", "\\2", input_results_folder)

  if(!tissue_code %in% tissue_db){
    stop("tissue_code: <", tissue_code, "> is not valid. Must be one of the followings:\n- ", paste(tissue_db, collapse = ", "))
  }else{
    return(tissue_code)
  }
}



