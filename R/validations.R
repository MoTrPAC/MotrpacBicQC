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
    processfolder <- stringr::str_extract(string = input_results_folder, pattern = "(BIC){0,1}RESULTS_[0-9]+")
  }

  if(is.na(processfolder)){
    stop("`PROCESS_YYYYMMDD` or `RESULTS_YYYYMMDD` folder is not recognize in the folder structure")
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
    stop("`BATCH#_YYYYMMDD` folder is not recognized in the folder structure.")
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



#' @title Extract PHASE from input folder path
#'
#' @description extract ASSAY from input folder path
#' @param input_results_folder (char) input_results_folder path
#' @param return_phase (char) return the phase only if `TRUE` (default)
#' @return (vector) PHASE code
#' @export
validate_phase <- function(input_results_folder, return_phase = TRUE){
  phase <- stringr::str_extract(string = input_results_folder,
                                pattern = "(PASS1A-06|PASS1A-18|PASS1B-06|PASS1B-18|PASS1C-06|PASS1C-18|HUMAN|HUMAN-PRECOVID|HUMAN-MAIN)")
  if(is.na(phase)){
    stop("Project phase (e.g. PASS1A-06) is not found in the folder structure, please, check guidelines")
  }else{
    if(return_phase) return(phase)
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



