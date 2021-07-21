# VALIDATIONS

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
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @export
check_viallabel_dmaqc <- function(vl_submitted,
                                  dmaqc_shipping_info,
                                  tissue_code,
                                  cas,
                                  phase,
                                  failed_samples,
                                  return_n_issues = FALSE,
                                  verbose = TRUE){
  
  # issue_count
  ic <- NA
  # There might be multiple phases to check: load both
  ph <- unlist(strsplit(phase, split = "\\|"))
  dmaqc_labels <- vector()
  for(i in 1:length(ph)){
    eph <- ph[i]
    pass <- gsub("(.*)(-)(.*)", "\\1", eph)
    month <- gsub("(.*)(-)(.*)", "\\3", eph)
    month <- as.integer(month)
    
    dmaqc_shipping_df <- read.delim(dmaqc_shipping_info, stringsAsFactors = FALSE)
    
    dmaqc_labels_temp <- dmaqc_shipping_df$vial_label[which(dmaqc_shipping_df$bic_tissue_code == tissue_code &
                                                              dmaqc_shipping_df$site_code == tolower(cas) &
                                                              dmaqc_shipping_df$phase == pass &
                                                              dmaqc_shipping_df$animal_age == month)]
    if(i == 1){
      dmaqc_labels <- as.character(dmaqc_labels_temp)
    }else{
      dmaqc_labels <- append(dmaqc_labels, as.character(dmaqc_labels_temp))
    }
  }
  
  if( length(dmaqc_labels) == 0){
    if(verbose) message("   + (+) DMAQC CHECK POINT: sample IDs not available in DMAQC dataset. Needs to be revised by BIC")
    ic <- "NOT_AVAILABLE"
  }else{
    if( setequal(vl_submitted, dmaqc_labels) ){
      if(verbose) message("   + (+) DMAQC CHECK POINT: samples sent to CAS have been processed: OK")
      ic <- "OK"
    }else{
      samples_missed <- setdiff(dmaqc_labels, vl_submitted)
      if(!is.null(failed_samples)){
        if(setequal(failed_samples, samples_missed)){
          if(verbose) message("   + (+) DMAQC CHECK POINT: samples sent to CAS have been processed (with known issues for some samples): OK")
          ic <- "OK"
        }else{
          if(verbose){
            message("      - (-) DMAQC CHECK POINT: samples not found in metadata_results: FAIL")
            message("\t - ", paste(samples_missed, collapse = "\n\t - "))
            ic <- "FAIL"
          }
        }
      }else{
        if(verbose){
          message("      - (-) DMAQC CHECK POINT: samples not found in metadata_results: FAIL")
          message("\t - ", paste(samples_missed, collapse = "\n\t - "))
        }
        ic <- "FAIL"
      }
    }
    
  }
  
  if(return_n_issues) return(ic)
}


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
                                       pattern = "(.*/BATCH\\d{1,2}\\_\\d{8})/")
  
  if(is.na(batch_folder)){
    stop("BATCH#_YYYYMMDD folder is not recognized in the folder structure.")
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
                                pattern = "(IONPNEG|RPNEG|RPPOS|HILICPOS|LRPPOS|LRPNEG|3HIB|AA|AC_DUKE|ACOA|BAIBA|CER_DUKE|CONV|KA|NUC|OA|SPHM|OXYLIPNEG|ETAMIDPOS|AC_MAYO|AMINES|CER_MAYO|TCA|PROT_PH|PROT_PR|PROT_AC|PROT_UB)")
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
                                pattern = "(PASS1A-06|PASS1A-18|PASS1B-06|PASS1B-18|HUMAN)")
  if(is.na(phase)){
    stop("<project phase (e.g. PASS1A-06)> is not found in the folder structure, please, check guidelines")
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
  tissue_code <- gsub("(.*)(T[0-9]{2,3})(.*)", "\\2", input_results_folder)

  if(!tissue_code %in% bic_animal_tissue_code$bic_tissue_code){
    stop("tissue_code: <", tissue_code, "> is not valid. Must be one of the following codes (check data object MotrpacBicQC::bic_animal_tissue_code):\n- ", paste(bic_animal_tissue_code$bic_tissue_code, collapse = ", "))
  }else{
    return(tissue_code)
  }
}



