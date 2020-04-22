

utils::globalVariables(
  c(".id",
    "bic_animal_tissue_code"))


# METABOLOMICS DATASETS: PRIMARY QC
#______________________________________________________________________________

#' @title check metadata metabolites
#'
#' @description check whether metadata_metabolites is following guidelines
#' @param df (data.frame) metadata_metabolites
#' @param nameun (char) specify whether `named` or `unnamed` files
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' check_metadata_metabolites(df = metadata_metabolites_named, nameun = "named")
#' }
#' @export
check_metadata_metabolites <- function(df,
                                       nameun,
                                       return_n_issues = FALSE,
                                       verbose = TRUE){

  # issue_count
  ic <- 0

  if(nameun == "named"){
    emeta_metabo_coln_named <- c("metabolite_name", "refmet_name", "rt", "mz", "neutral_mass", "formula")
  }else if(nameun == "unnamed"){
    emeta_metabo_coln_named <- c("metabolite_name", "rt", "mz", "neutral_mass")
  }else{
    stop("{nameun} option not valid. Options: named/unnamed")
  }

  colnames(df) <- tolower(colnames(df))

  if(all(emeta_metabo_coln_named %in% colnames(df))){
    if(verbose) message("   + (+) All required columns present")
  }

  # Evaluate every column
  flag_mm <- FALSE
  if("metabolite_name" %in% colnames(df)){
    if(length(unique(df$metabolite_name)) != dim(df)[1]){
      duplis_details <- df$metabolite_name[duplicated(df$metabolite_name)]
      duplis <- length(duplis_details)
      if(verbose) message("      - (-) {metabolite_name} non-unique values detected: ", duplis)
      if(verbose) message("\n\t\t - ", paste(duplis_details, collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      flag_mm <- TRUE
      if(verbose) message("   + (+) {metabolite_name} OK")
    }
  }else{
    if(verbose) message("      - (-) {metabolite_name} is missed: FAIL")
    ic <- ic + 1
  }

  # Metabolomics data dictionary
  mdd <- get_and_validate_mdd()

  # refmet_name only expected on named metabolites
  if(nameun == "named"){
    if("refmet_name" %in% colnames(df)){
      if(length(unique(df$refmet_name)) != dim(df)[1]){
        duplis_details <- df$refmet_name[duplicated(df$refmet_name)]
        duplis <- length(duplis_details)
        if(verbose) message("      - (-) {refmet_name} non-unique values detected: ", duplis)
        if(verbose) message("\n\t\t - ", paste(duplis_details, collapse = "\n\t\t - "))
        ic <- ic + 1
      }else{
        if(verbose) message("   + (+) {refmet_name} unique values: OK")
      }

      if(!all(df$refmet_name %in% mdd$refmet_name)){
        if(verbose) message("      - (-) {refmet_name}: ids not found in RefMet Metabolomics Data Dictionary: FAIL")
        if(verbose) message("\n\t\t - ", paste(setdiff(df$refmet_name, mdd$refmet_name), collapse = "\n\t\t - "))
        ic <- ic + 1
      }else{
        if(verbose) message("   + (+) {refmet_name} ids found in refmet: OK")
      }

    }else{
      if(verbose) message("      - (-) {refmet_name} column missed: FAIL")
      ic <- ic + 1
    }
  }

  if("rt" %in% colnames(df)){
    if(!all(is.numeric(df$rt))){
      if(verbose) message("      - (-) {rt} non numeric values found: FAIL")
      nonnum <- df$rt[which(!grepl('^[0-9]', df$rt))]
      if(verbose) message("\n\t\t - ", paste(nonnum, collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("   + (+) {rt} all numeric: OK")
    }
  }else{
    if(verbose) message("      - (-) {rt} column missed: FAIL")
    ic <- ic + 1
  }

  if("mz" %in% colnames(df)){
    if(!all(is.numeric(df$mz))){
      if(verbose) message("      - (-) {mz}: non numeric values found: FAIL")
      nonnum <- df$mz[which(!grepl('^[0-9]', df$mz))]
      if(verbose) message("\n\t\t - ", paste(nonnum, collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("   + (+) {mz} all numeric: OK")
    }
  }else{
    if(verbose) message("      - (-) {mz} column missed: FAIL")
    ic <- ic + 1
  }

  if(nameun == "named"){
    if("neutral_mass" %in% colnames(df)){
      if(!all(is.numeric(df$neutral_mass))){
        if(verbose) message("      - (-) {neutral_mass}: non numeric values found: FAIL")
        nonnum <- df$neutral_mass[which(!grepl('^[0-9]', df$neutral_mass))]
        if(verbose) message("\n\t\t - ", paste(nonnum, collapse = "\n\t\t - "))
        ic <- ic + 1
      }else{
        if(verbose) message("   + (+) {neutral_mass} all numeric values OK")
      }
    }else{
      if(verbose) message("      - (-) {neutral_mass} column missed: FAIL")
      ic <- ic + 1
    }
    # Evaluate every column
    if(!("formula" %in% colnames(df))){
      if(verbose) message("      - (-) {formula} column is missed: FAIL")
      ic <- ic + 1
    }else{
      if(verbose) message("   + (+) {formula} available: OK")
    }
  }else{
    if("neutral_mass" %in% colnames(df)){
      if(!all(is.numeric(df$neutral_mass))){
        if(verbose) message("      - (-) {neutral_mass}: non numeric values found: FAIL")
        nonnum <- df$neutral_mass[which(!grepl('^[0-9]', df$neutral_mass))]
        if(verbose) message("\n\t\t - ", paste(nonnum, collapse = "\n\t\t - "))
        ic <- ic + 1
      }else{
        if(verbose) message("   + (+) {neutral_mass} all numeric values OK")
      }
    }else{
      if(verbose) message("   + ( ) {neutral_mass} column missed (but not-required for UNNAMED)")
    }
  }

  if(return_n_issues) return(ic)

} #end check_metadata_metabolites



# ------------------------------------------------------------------------------
#' @title check metadata samples
#'
#' @description check whether metadata_sample is following guidelines
#' @param df (data.frame) metadata_metabolites
#' @param cas (char) CAS site code
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' check_metadata_samples(df = metadata_sample_named, cas = "umichigan")
#' }
#' @export
check_metadata_samples <- function(df,
                                   cas,
                                   return_n_issues = FALSE,
                                   verbose = TRUE){

  # issue_count
  ic = 0

  # Columns to check
  emeta_sample_coln <- c("sample_id", "sample_type", "sample_order", "raw_file")

  colnames(df) <- tolower(colnames(df))

  if( all(emeta_sample_coln %in% colnames(df)) ){
    if(verbose) message("   + (+) All required columns present")
  }else{
    if(verbose) message("      - (-) Expected COLUMN NAMES are missed: FAIL")
    ic <- ic + 1
  }

  # Check every column

  # sample_id: si
  flag_si <- FALSE
  if( "sample_id" %in% colnames(df) ){
    if( length(unique(df$sample_id)) != dim(df)[1] ){
      if(verbose) message("      - (-) {sample_id}: Non-unique values detected: FAIL")
      ic <- ic + 1
    }else{
      flag_si <- TRUE
      if(verbose) message("   + (+) {sample_id} seems OK")
    }
  }else{
    if(verbose) message("      - (-) {metabolite_name} is missed: FAIL")
    ic <- ic + 1
  }

  # sample_type: st
  esample_types <- c("Sample", "QC-Pooled", "QC-Reference", "QC-Blank",
                     "QC-Identification", "QC-InternalStandard", "QC-PreRun",
                     "QC-ExternalStandard", "QC-DriftCorrection")

  if("sample_type" %in% colnames(df)){
    if(!all(df$sample_type %in% esample_types)){
      if(verbose) message("      - (-) Error: undefined sample types: ", appendLF = FALSE)
      if(verbose) message("\n\t\t - ", paste(setdiff(df$sample_type, esample_types), collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("   + (+) {sample_type} seems OK")
    }
  }else{
    if(verbose) message("      - (-) {refmet_name} column missed: FAIL")
    ic <- ic + 1
  }


  if("sample_order" %in% colnames(df)){
    if(!all(is.numeric(df$sample_order))){
      if(verbose) message("      - (-) {sample_order} non numeric values found -> ", appendLF = FALSE)
      nonnum <- df$sample_order[which(!grepl('^[0-9]', df$sample_order))]
      if(verbose) message("\n\t\t - ", paste(nonnum, collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("   + (+) {sample_order} is numeric")
    }

    if(length(unique(df$sample_order)) != dim(df)[1]){
      if(verbose) message("      - (-) {sample_order}: Non-unique values detected -> ", appendLF = FALSE)
      if(verbose) message(paste(df$sample_order[unique(duplicated(df$sample_order))]))
      ic <- ic + 1
    }else{
      if(verbose) message("   + (+) {sample_order} unique values OK")
    }

  }else{
    if(verbose) message("      - (-) {sample_order} column missed: FAIL")
    ic <- ic + 1
  }

  if("raw_file" %in% colnames(df)){
    if( cas != "emory"){
      if(length(unique(df$raw_file)) != dim(df)[1]){
        if(verbose) message("      - (-) {raw_file}: Non-unique values detected or missed -> ", appendLF = FALSE)
        if(verbose) message(paste( unique(df$raw_file[unique(duplicated(df$raw_file))]) ) )
        ic <- ic + 1
      }else{
        if(verbose) message("   + (+) {raw_file} unique values OK")
      }
    }
  }else{
    if(verbose) message("      - (-) {raw_file} column missed: FAIL")
    ic <- ic + 1
  }

  if(return_n_issues) return(ic)
} #end check_metadata_samples


# RESULTS

# ------------------------------------------------------------------------------
#' @title check results
#'
#' @description check whether results file is following guidelines
#' @param r_m (data.frame) results df
#' @param m_s (char) metadata_sample df
#' @param m_m (char) metadata_metabolites df
#' @param return_n_issues (logical) if `TRUE` returns the number of issues
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' check_results(r_m = results_named, m_s = metadata_sample_named, m_m = metadata_metabolites_named)
#' }
#' @export
check_results <- function(r_m,
                          m_s,
                          m_m,
                          return_n_issues = FALSE,
                          verbose = TRUE){

  # issue_count
  ic = 0

  # Expected Columns that should be present
  eresults_coln <- c("metabolite_name", unique(m_s$sample_id))

  flag_out <- TRUE
  if(!setequal(colnames(r_m), eresults_coln)){
    extra_in_results <- setdiff(colnames(r_m), eresults_coln)
    if(length(extra_in_results > 0)){
      if(verbose) message("\n      - (-) Column(s) NOT expected in {results_metabolite} file which are missed in {metadata_samples}: \n\t\t - ", paste(extra_in_results, collapse = "\n\t\t - "))
    }

    extra_in_msr <- setdiff(eresults_coln, colnames(r_m))
    if(length(extra_in_msr)){
      if(verbose) message("\n      - (-) Column(s) available in {metadata_samples} missed in {results_metabolite}: \n\t\t - ", paste(extra_in_msr, collapse = "\n\t\t - "))
    }
    flag_out <- FALSE
    ic <- ic + 1
  }else{
    if(verbose) message("   + (+) All samples from [results_metabolite] are available in [metadata_sample]")
  }

  # Check that metabolites names matches
  if(!is.null(m_m)){
    if("metabolite_name" %in% colnames(r_m) & "metabolite_name" %in% colnames(m_m)){
      if(!setequal(r_m$metabolite_name, m_m$metabolite_name)){
        if(verbose) message("      - (-) {metabolite_name} in [results], not found in [metadata_metabolites]:\n\t\t- ",
                paste(setdiff(r_m$me, m_m$metabolite_name), collapse = "\n\t\t- "))
        ic <- ic + 1
      }else{
        if(verbose) message("   + (+) {metabolite_name} is identical in both [results] and [metadata_metabolites] files: OK")
      }
    }
  }

  # Check if sample columns are numeric (but only if sample matches between them)
  if(flag_out){
    r_m_num <- r_m[,m_s$sample_id] %>% dplyr::select_if(is.numeric)
    if( !identical(r_m_num, r_m[,m_s$sample_id]) ){
      if(verbose) message("      - (-) Non-numeric columns identified")
      r_m_nn <- r_m[,m_s$sample_id] %>% dplyr::select_if(negate(is.numeric))
      if(verbose) message("      - (-) Non-numeric columns: ", paste(colnames(r_m_nn), collapse = ","))
      ic <- ic + 1
    }else{
      if(verbose) message("   + (+) {sample_id} columns are numeric: OK")
    }
  }else{
    ic <- ic + 1
  }

  if(return_n_issues) return(ic)
} # check_results

# ------------------------------------------------------------------------------
#' @title check rawfiles between manifest and metabolite_sample matches
#'
#' @description check rawfiles between manifest and metabolite_sample matches
#' @param input_results_folder (char) input path folder
#' @param m_s_n_raw (list) list of raw files available in the metadata sample file
#' @param return_n_issues (logical) if `TRUE` returns the number of issues
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @export
check_manifest_rawdata <- function(input_results_folder,
                                   m_s_n_raw,
                                   return_n_issues = FALSE,
                                   verbose = TRUE){

  # issue_count
  ic <- 0

  # get the raw file location
  input_results_raw <- gsub("PROCESSED_\\d+", "RAW", input_results_folder)
  filepattern <- "manifest"

  # Get file matching pattern
  file_rawmanifest <- list.files(input_results_raw,
                                 pattern=filepattern,
                                 full.names=TRUE,
                                 recursive = TRUE,
                                 ignore.case = TRUE)

  # Check if file is found and deal with many files
  if(length(file_rawmanifest) != 1){
    if(length(file_rawmanifest) >= 1){
      if(verbose) message("      - (-) open_file: more than one file detected")
      if(verbose) message("\n\t\t - ", paste(file_rawmanifest, collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("   + (-) File [", filepattern, "] not found")
    }
    flag <- FALSE
    ofile <- NULL
  }else{
    flag <- TRUE
    ofile <- read.delim(file_rawmanifest[1], stringsAsFactors = FALSE, check.names = FALSE)
  }

  if(flag){
    if(nrow(ofile) == 0){
      if(verbose) message("   + ( ) MANIFEST raw file is empty")
      return(NULL)
    }else{
      if("raw_file" %in% colnames(ofile)){
        if(setequal(ofile$raw_file, m_s_n_raw)){
          if(verbose) message ("   + (+) RAW files match between both [metadata_samples] and [raw/manifest] files: OK")
        }else{
          if(verbose) message("      - ( ) RAW files DO NOT match between both [metadata_samples] and [raw/manifest] files. Investigating:")

          extra_in_a <- setdiff(ofile$raw_file, m_s_n_raw)
          if(length(extra_in_a) > 0){
            if(verbose) message("      - ", length(extra_in_a)," extra {raw_file} found in [raw/manifest]. This is: OK")
          }

          extra_in_b <- setdiff(m_s_n_raw, ofile$raw_file)
          if(length(extra_in_b) > 0){
            if(verbose) message("      - (-) {raw_file}(s) only found in [metadata_samples] not available in [manifest] file: FAIL\n\t\t - ", paste(extra_in_b, collapse = "\n\t\t - "))
            ic <- ic + 1
          }
        }
      }else{
        if(verbose) message("      - (-) {raw_file} column not found")
        ic <- ic + 1
      }
      if("md5" %in% colnames(ofile)){
        if(verbose) message("   + (+) {md5} column found: OK")
      }else{
        if(verbose) message("      - (-) {md5} column not found: FAIL")
        ic <- ic + 1
      }
    }
  }
  if(return_n_issues) return(ic)
}

# ------------------------------------------------------------------------------
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
      if(verbose) message("      - (-) open_file: more than one file detected: FAIL")
      if(verbose) message("\n\t\t - ", paste(file_metametabolites, collapse = "\n\t\t - "))
    }else{
      if(verbose) message("   + ( ) File [", filepattern, "] not found")
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
        return(ofile$sample_id)
      }else{
        if(verbose) message("      - (-) {sample_id} column not found: FAIL")
        return(NULL)
      }
    }
  }
}

# ------------------------------------------------------------------------------
#' @title validate vial label from DMAQC
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

  pass <- gsub("(.*)(-)(.*)", "\\1", phase)
  month <- gsub("(.*)(-)(.*)", "\\3", phase)

  dmaqc_labels <- dmaqc_shipping_info$vial_label[which(dmaqc_shipping_info$bic_tissue_code == tissue_code &
                                                         dmaqc_shipping_info$site_code == tolower(cas) &
                                                         dmaqc_shipping_info$phase == pass &
                                                         dmaqc_shipping_info$animal_age == month)]
  if(setequal(vl_submitted, dmaqc_labels)){
    if(verbose) message("   + (+) DMAQC CHECK POINT: samples sent to CAS have been processed: OK")
    ic <- "OK"
  }else{
    samples_missed <- setdiff(dmaqc_labels, vl_submitted)
    if(!is.null(failed_samples)){
      if(setequal(failed_samples, samples_missed)){
        if(verbose) message("   + (+) DMAQC CHECK POINT: samples sent to CAS have been processed (with known issues for some samples): OK")
        ic <- "OK"
      }
    }else{
      if(verbose) message("      - (-) DMAQC CHECK POINT: samples not found in metadata_results: FAIL")
      if(verbose) message("\t - ", paste(samples_missed, collapse = "\n\t - "))
      ic <- "FAIL"
    }

  }
  if(return_n_issues) return(ic)
}

# ------------------------------------------------------------------------------
#' @title Validate a Metabolomics submission
#'
#' @description Validate a Metabolomics submission
#' @param input_results_folder (char) path to the PROCESSED folder to check
#' @param tissue_code (char) tissue code
#' @param cas (char) CAS code
#' @param phase (char) phase code
#' @param dmaqc_shipping_info (char) phase code
#' @param return_n_issues (logical) if `TRUE` returns the number of issues
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (data.frame) Summary of issues
#' @export
validate_metabolomics <- function(input_results_folder,
                                  cas,
                                  tissue_code,
                                  phase,
                                  dmaqc_shipping_info,
                                  return_n_issues = FALSE,
                                  verbose = TRUE){

  # issue_count
  ic <- 0
  ic_m_m_n <- NA
  ic_m_m_u <- NA
  ic_m_s_n <- NA
  ic_m_s_u <- NA
  ic_r_m_n <- NA
  ic_r_m_u <- NA
  ic_mrd <- NA
  ic_vl <- NA
  vial_label <- NA
  qc_samples <- NA

  input_folder_short <- regmatches(input_results_folder, regexpr("PASS.*PROCESSED_[0-9]{8}", input_results_folder))
  if(is_empty(input_folder_short)){
    if(verbose) message("\nThe PROCESSED_YYYYMMDD folder full path is not correct. Example:")
    if(verbose) message("/full/path/to/folder/PASS1A-06/T66/RPNEG/BATCH1_20190822/PROCESSED_202003")
    stop("Input folder not according to guidelines")
  }

  if(verbose) message("# METABOLOMICS QC report\n\n")
  if(verbose) message("+ Site: ", cas)
  if(verbose) message("+ Folder: `",paste0(input_folder_short),"`")

  # Is a targeted site? Unname compounds not checked
  untargeted <- TRUE

  if(cas %in% c("mayo", "emory", "duke")){
    untargeted <- FALSE
  }

  if(verbose) message("\n## QC metadata_metabolites")
  cat("\n")
  if(verbose) message("*NAMED metadata metabolites*\n")

  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata_metabolites_named.*.txt|Metadata_named_metabolites.txt")
  f_mmn <- lista$flag
  if(f_mmn){
    m_m_n <- lista$df
    ic_m_m_n <- check_metadata_metabolites(m_m_n, "named", return_n_issues = TRUE)
  }else{
    if(verbose) message("      - (-) {metadata_metabolites_named} not available")
    ic <- ic + 1
  }

  if(untargeted){
    if(verbose) message("\n*UNNAMED metadata metabolites*\n")
    lista <- open_file(input_results_folder, "metadata_metabolites_unnamed.*.txt|Metadata_unnamed_metabolites.txt")
    f_mmu <- lista$flag
    if(f_mmu) {
      m_m_u <- lista$df
      ic_m_m_u <- check_metadata_metabolites(m_m_u, "unnamed", return_n_issues = TRUE)
    }else{
      if(verbose) message("      - (-) {metadata_metabolites_unnamed} not available")
      ic <- ic + 1
    }
  }

  if(verbose) message("\n\n## QC metadata_sample files\n")

  if(verbose) message("\n*NAMED metadata_sample*\n")

  lista <- open_file(input_results_folder, filepattern = "metadata_sample.*_named.*.txt")
  f_msn <- lista$flag
  if(f_msn){
    m_s_n <- lista$df
    ic_m_s_n <- check_metadata_samples(df = m_s_n, cas = cas, return_n_issues = TRUE)
    # Extract the number of samples
    if(!is.null(m_s_n)){
      #Double check that the columns are there
      if( all(c("sample_id", "sample_type") %in% colnames(m_s_n)) ){
        vial_label <- length(m_s_n$sample_id[which(m_s_n$sample_type == "Sample")])
        qc_samples <- length(m_s_n$sample_id[which(m_s_n$sample_type != "Sample")])
      }
    }
  }

  if(untargeted){
    if(verbose) message("\n*UNNAMED metadata_sample*\n")
    lista <- open_file(input_results_folder, "metadata_sample.*_unnamed.*.txt")
    f_msu <- lista$flag
    if(f_msu){
      m_s_u <- lista$df
      ic_m_s_u <- check_metadata_samples(m_s_u, cas, return_n_issues = TRUE)
    }

    # NAMED AND UNNAMED MUST MATCH TESTS
    if(f_msn & f_msu){
      if(isTRUE(all.equal(m_s_n, m_s_u))){
        if(verbose) message("   + (+) Metadata samples: named and unnamed are identical: OK")
      }else{
        if(verbose) message("      - (-) Metadata samples: named and unnamed files differ")
        ic <- ic + 1
      }
    }
  }

  if(verbose) message("\n\n## QC Results\n")

  if(verbose) message("\n*NAMED results_metabolites*\n")

  # if(cas == "emory"){
  #   lista <- open_file(input_results_folder, "results_metabolites_named.*adjusted.*.txt")
  # }else{
  lista <- open_file(input_results_folder, "results_metabolites_named.*.txt|Results_named_metabolites.txt")
  # }

  f_rmn <- lista$flag
  if(f_rmn){
    r_m_n <- lista$df
    if(f_msn & f_mmn){
      ic_r_m_n <- check_results(r_m = r_m_n,
                                m_s = m_s_n,
                                m_m = m_m_n,
                                return_n_issues = TRUE)
    }
  }else{
    if(verbose) message("      - (-) RESULTS-NAMED file not available")
    ic <- ic + 1
  }

  if(untargeted){
    if(verbose) message("\n*UNNAMED results_metabolites*\n")
    lista <- open_file(input_results_folder, "results_metabolites_unnamed.*.txt|Results_unnamed_metabolites.txt")
    f_rmu <- lista$flag
    if(f_rmu){
      r_m_u <- lista$df
      if(f_msu & f_mmu){
        ic_r_m_u <- check_results(r_m_u,
                                  m_s_u,
                                  m_m_u,
                                  return_n_issues = TRUE)
      }
    }else{
      if(verbose) message("      - (-) UNNAMED-RESULTS file not available ")
      ic <- ic + 1
    }

    if(f_rmn & f_rmu){
      if(verbose) message("- Checking columns match between named and unnamed results\n")
      if(setequal(colnames(r_m_n), colnames(r_m_u))){
        if(verbose) message("   + (+) All samples found on both NAMED and UNNAMED files: OK")
      }else{
        if(verbose) message("\n      - (-) Samples DO NOT MATCH on both NAMED and UNNAMED files")
        ic <- ic + 1
        extra_in_rnamed <- setdiff(colnames(r_m_n), colnames(r_m_u))
        if(length(extra_in_rnamed) > 0){
          if(verbose) message("\n      - (-) Column(s) only found in [results NAMED] file:\n", paste(extra_in_rnamed, collapse = "\n\t\t - "))
        }

        extra_in_unamed <- setdiff(colnames(r_m_u), colnames(r_m_n))
        if(length(extra_in_unamed) > 0){
          if(verbose) message("\n      - (-) Column(s) only found in [results UNNAMED] file:\n", paste(extra_in_unamed, collapse = "\n\t\t - "))
        }
      }
    }
  }

  if(f_msn){
    if(verbose) message("\n\n## QC raw_file information\n")
    ic_mrd <- check_manifest_rawdata(input_results_folder = input_results_folder,
                                     m_s_n_raw = unique(m_s_n$raw_file),
                                     return_n_issues = TRUE)
  }else{
    if(verbose) message("\n\n## Validate {raw_files} match between [RAW/Manifest.txt] and [metadata_samples]\n")
    if(verbose) message("      (-) FAIL (medatada_samples file not available)")
    ic <- ic + 1
  }

  if(verbose) message("\n\n## DMAQC validation\n")
  failed_samples <- check_failedsamples(input_results_folder = input_results_folder)

  # Validate vial labels
  if(f_msn){
    vl_results <- m_s_n$sample_id[which(m_s_n$sample_type == "Sample")]
    ic_vl <- check_viallabel_dmaqc(vl_submitted = vl_results,
                                   tissue_code = tissue_code,
                                   cas = cas,
                                   phase = phase,
                                   failed_samples = failed_samples,
                                   dmaqc_shipping_info = dmaqc_shipping_info,
                                   return_n_issues = TRUE)
  }

  batchversion <- stringr::str_extract(string = input_results_folder, pattern = "BATCH.*_[0-9]+/PROCESSED_[0-9]+")
  assay <- stringr::str_extract(string = input_results_folder,
                                pattern = "IONPNEG|RPNEG|RPPOS|HILICPOS|LRPPOS|LRPNEG|OXYLIPNEG")

  qc_date <- Sys.time()
  qc_date <- gsub("-", "", qc_date)
  qc_date <- gsub(" ", "_", qc_date)
  qc_date <- gsub(":", "", qc_date)

  reports <- data.frame(cas = cas,
                        phase= phase,
                        tissue = tissue_code,
                        t_name = bic_animal_tissue_code$bic_tissue_name[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)],
                        assay = assay,
                        version = batchversion,
                        vial_label = vial_label,
                        qc_samples = qc_samples,
                        dmaqc_valid = ic_vl,
                        critical_issues = ic,
                        raw_manifest = ic_mrd,
                        m_metab_n = ic_m_m_n,
                        m_metab_u = ic_m_m_u,
                        m_sample_n = ic_m_s_n,
                        m_sample_u = ic_m_s_u,
                        results_n = ic_r_m_n,
                        results_u = ic_r_m_u,
                        qc_date = qc_date)

  if(return_n_issues) return(reports)
}





