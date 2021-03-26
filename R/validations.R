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
    processfolder <- stringr::str_extract(string = input_results_folder, pattern = "RESULTS_[0-9]+")
  }

  if(is.na(processfolder)){
    stop("PROCESS_YYYYMMDD or RESULTS_YYYYMMDD folder is not recognize in the folder structure")
  }else{
    return(processfolder)
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
                                       pattern = "(.*/BATCH\\d{1}\\_\\d{8})/")
  
  if(is.na(batch_folder)){
    stop("BATCH#_YYYYMMDD folder is not recognized in the folder structure")
  }else{
    return(batch_folder)
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
                                pattern = "(IONPNEG|RPNEG|RPPOS|HILICPOS|LRPPOS|LRPNEG|3HIB|AA|AC_DUKE|ACOA|BAIBA|CER_DUKE|CON|KA|NUC|OA|SPHM|OXYLIPNEG|ETAMIDPOS|AC_MAYO|AMINES|CER_MAYO|TCA|PROT_PH|PROT_PR|PROT_AC|PROT_UB)")
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
  tissue_code <- gsub("(.*)(T[0-9]{2})(.*)", "\\2", input_results_folder)

  if(!tissue_code %in% bic_animal_tissue_code$bic_tissue_code){
    stop("tissue_code: <", tissue_code, "> is not valid. Must be one of the following codes (check data object MotrpacBicQC::bic_animal_tissue_code):\n- ", paste(bic_animal_tissue_code$bic_tissue_code, collapse = ", "))
  }else{
    return(tissue_code)
  }
}



