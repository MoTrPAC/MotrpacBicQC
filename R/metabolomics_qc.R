

utils::globalVariables(
  c("bic_animal_tissue_code",
    "phenotypes_pass1a06_short"))


# METABOLOMICS DATASETS: PRIMARY QC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title check metadata metabolites
#'
#' @description check whether metadata_metabolites is following guidelines
#' @param df (data.frame) metadata_metabolites
#' @param name_id (char) specify whether `named` or `unnamed` files
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' check_metadata_metabolites(df = metadata_metabolites_named, name_id = "named")
#' }
#' @export
check_metadata_metabolites <- function(df,
                                       name_id,
                                       return_n_issues = FALSE,
                                       verbose = TRUE){

  # issue_count
  ic <- 0

  if(name_id == "named"){
    emeta_metabo_coln_named <- c("metabolite_name", "refmet_name", "rt", "mz", "neutral_mass", "formula")
  }else if(name_id == "unnamed"){
    emeta_metabo_coln_named <- c("metabolite_name", "rt", "mz", "neutral_mass")
  }else{
    stop("{name_id} option not valid. Options: named/unnamed")
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
  if(name_id == "named"){
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

  if(name_id == "named"){
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

  # colnames(df) <- tolower(colnames(df))

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Validate a Metabolomics submission
#'
#' @description Validate a Metabolomics submission
#' @param input_results_folder (char) path to the PROCESSED folder to check
#' @param cas (char) CAS code
#' @param dmaqc_shipping_info (char) phase code
#' @param return_n_issues (logical) if `TRUE` returns the number of issues
#' @param full_report (logical) if `FALSE` (default) it returns only the
#' total number of issues. If `TRUE` returns the details of the number of issues (by
#' group of files, e.g., results, metadata_metabolites, etc)
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (data.frame) Summary of issues
#' @export
validate_metabolomics <- function(input_results_folder,
                                  cas,
                                  dmaqc_shipping_info = NULL,
                                  return_n_issues = FALSE,
                                  full_report = FALSE,
                                  verbose = TRUE){

  validate_cas(cas = cas)

  processfolder <- validate_processFolder(input_results_folder)

  assay <- validate_assay(input_results_folder)

  phase <- validate_phase(input_results_folder)

  tissue_code <- validate_tissue(input_results_folder)

  # issue_count
  ic <- 0
  ic_m_m_n <- NA
  ic_m_m_u <- NA
  ic_m_s_n <- NA
  ic_m_s_u <- NA
  ic_r_m_n <- NA
  ic_r_m_u <- NA
  ic_mrd <- NA
  if(is.null(dmaqc_shipping_info)){
    ic_vl <- "missed"
  } else{
    ic_vl <- NA
  }
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
  r_m_n <- NULL
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

  # Validate vial labels from DMAQC
  if(is.na(ic_vl)){
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
  }

  if(ic > 4){
    message("\nERROR: Too many files missed. Please, provide the full path to PROCESSED_YYYYMMDD folder
         input_results_folder = /full/path/to/PROCESSED_YYYYMMDD\n")
    stop("Revise <input_results_folder> argument")
  }

  batchversion <- stringr::str_extract(string = input_results_folder, pattern = "BATCH.*_[0-9]+/PROCESSED_[0-9]+")

  qc_date <- Sys.time()
  qc_date <- gsub("-", "", qc_date)
  qc_date <- gsub(" ", "_", qc_date)
  qc_date <- gsub(":", "", qc_date)
  t_name <- bic_animal_tissue_code$bic_tissue_name[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)]

  if(return_n_issues){
    issues <- sum(ic, ic_mrd, ic_m_m_n, ic_m_m_u, ic_m_s_n, ic_m_s_u, ic_r_m_n, ic_r_m_u, na.rm = TRUE)
    if(verbose) message("\nTOTAL NUMBER OF ISSUES: ", issues,"\n")
    if(full_report){
      reports <- data.frame(cas = cas,
                            phase= phase,
                            tissue = tissue_code,
                            t_name = t_name,
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
      return(reports)
    }
  }
}

#' @title Combines all files from a Metabolomics batch
#'
#' @description Combines all the files from an untargeted assay submission
#' including metadata_sample, metadata_metabolites, results for both "NAMED" and "UNNAMED"
#' folders
#' @param input_results_folder (char) Path to the PROCESSED_YYYYMMDD folder
#' @param cas (char) Chemical Analytical Site code (e.g "umichigan")
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples
#' \dontrun{
#' all_datasets <- combine_metabolomics_batch(
#'                         input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD/",
#'                         cas = "umichigan")
#' }
#' @export
combine_metabolomics_batch <- function(input_results_folder,
                                       cas,
                                       verbose = TRUE){

  m_m_n = m_m_u = m_s_n = r_m_n = r_m_u = NULL

  # Validations----
  cas <- tolower(cas)
  validate_cas(cas)
  phase <- validate_phase(input_results_folder)

  if(grepl("^PASS1A-06$", phase)){
    stop("<combine_metabolomics_batch> only supports PASS1A batches")
  }

  processfolder <- validate_processFolder(input_results_folder)
  assay <- validate_assay(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)

  # Output name----
  output_name <- paste0(cas, ".", tissue_code, ".", tolower(phase), ".",tolower(assay), ".", tolower(processfolder))

  # COUNTS AND FLAGS
  ic <- 0

  vial_label <- NA
  qc_samples <- NA

  # Load Metabolomics
  if(verbose) message("# LOAD METABOLOMICS BATCH")
  if(verbose) message("+ Site: ", cas)
  if(verbose) message("+ Folder: `",paste0(input_results_folder),"`")

  # Is a targeted site? Unname compounds not checked
  untargeted <- TRUE

  if(cas %in% c("mayo", "emory", "duke")){
    untargeted <- FALSE
  }

  # metadata_metabolites------
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata_metabolites_named.*.txt|Metadata_named_metabolites.txt",
                     verbose = verbose)
  f_mmn <- lista$flag
  if(f_mmn){
    m_m_n <- lista$df
    check_metadata_metabolites(m_m_n, "named",
                               return_n_issues = FALSE,
                               verbose = verbose)
  }else{
    if(verbose) message("      - (-) {metadata_metabolites_named} not available")
    ic <- ic + 1
  }

  if(untargeted){
    if(verbose) message("\n*UNNAMED metadata metabolites*\n")
    lista <- open_file(input_results_folder,
                       "metadata_metabolites_unnamed.*.txt|Metadata_unnamed_metabolites.txt",
                       verbose = verbose)
    f_mmu <- lista$flag
    if(f_mmu) {
      m_m_u <- lista$df
      check_metadata_metabolites(m_m_u, "unnamed",
                                 return_n_issues = FALSE,
                                 verbose = verbose)
    }else{
      if(verbose) message("      - (-) {metadata_metabolites_unnamed} not available")
      ic <- ic + 1
    }
  }

  # metadata_sample-----
  if(verbose) message("\n\n## Metadata_sample files\n")
  if(verbose) message("\n*NAMED metadata_sample*\n")

  lista <- open_file(input_results_folder,
                     "metadata_sample_named.*.txt|Metadata_named_samples.txt",
                     verbose = verbose)
  f_msn <- lista$flag
  if(f_msn){
    m_s_n <- lista$df
    check_metadata_samples(df = m_s_n,
                           cas = cas,
                           return_n_issues = FALSE,
                           verbose = verbose)
    # Extract the number of samples
    if(!is.null(m_s_n)){
      #Double check that the columns are there
      if( all(c("sample_id", "sample_type") %in% colnames(m_s_n)) ){
        vial_label <- length(m_s_n$sample_id[which(m_s_n$sample_type == "Sample")])
        qc_samples <- length(m_s_n$sample_id[which(m_s_n$sample_type != "Sample")])
        if(verbose) message("\t+ Number of Samples: ",vial_label, " QCs: ", qc_samples)
      }
    }
  }

  if(untargeted){
    if(verbose) message("\n*UNNAMED metadata_sample*\n")
    lista <- open_file(input_results_folder,
                       "metadata_sample_unnamed.*.txt|Metadata_unnamed_samples.txt",
                       verbose = verbose)
    f_msu <- lista$flag
    if(f_msu){
      m_s_u <- lista$df
      check_metadata_samples(m_s_u, cas,
                             return_n_issues = FALSE,
                             verbose = verbose)
    }

    # NAMED AND UNNAMED MUST MATCH TESTS
    if(f_msn & f_msu){
      if(isTRUE(all.equal(m_s_n, m_s_u))){
        if(verbose) message("   + (+) Metadata samples: named and unnamed are identical: OK")
      }else{
        stop("      - (-) Metadata samples: named and unnamed files differ")
        # ic <- ic + 1
      }
    }
  }


  # results --------------------------------------------------------------------
  if(verbose) message("\n\n## Results\n")
  if(verbose) message("\n*NAMED results_metabolites*\n")

  lista <- open_file(input_results_folder,
                     "results_metabolites_named.*.txt|Results_named_metabolites.txt",
                     verbose = verbose)

  f_rmn <- lista$flag
  if(f_rmn){
    r_m_n <- lista$df
    if(f_msn & f_mmn){
      check_results(r_m = r_m_n,
                    m_s = m_s_n,
                    m_m = m_m_n,
                    return_n_issues = FALSE,
                    verbose = verbose)
    }
  }else{
    if(verbose) message("      - (-) RESULTS-NAMED file not available")
    ic <- ic + 1
  }

  if(untargeted){
    if(verbose) message("\n*UNNAMED results_metabolites*\n")
    lista <- open_file(input_results_folder,
                       "results_metabolites_unnamed.*.txt|Results_unnamed_metabolites.txt",
                       verbose = verbose)
    f_rmu <- lista$flag
    if(f_rmu){
      r_m_u <- lista$df
      if(f_msu & f_mmu){
        check_results(r_m_u,
                      m_s_u,
                      m_m_u,
                      return_n_issues = FALSE,
                      verbose = verbose)
      }
    }else{
      if(verbose) message("      - (-) UNNAMED-RESULTS file not available ")
      ic <- ic + 1
    }

    if(f_rmn & f_rmu){
      if(verbose) message("- Checking columns match between named and unnamed results")
      if(setequal(colnames(r_m_n), colnames(r_m_u))){
        if(verbose) message("   + (+) All samples found on both NAMED and UNNAMED files: OK")
      }else{
        if(verbose) message("      - (-) Samples DO NOT MATCH on both NAMED and UNNAMED files")
        ic <- ic + 1
        extra_in_rnamed <- setdiff(colnames(r_m_n), colnames(r_m_u))
        if(length(extra_in_rnamed) > 0){
          if(verbose) message("      - (-) Column(s) only found in [results NAMED] file: ", paste(extra_in_rnamed, collapse = "\n\t\t - "))
        }

        extra_in_unamed <- setdiff(colnames(r_m_u), colnames(r_m_n))
        if(length(extra_in_unamed) > 0){
          if(verbose) message("      - (-) Column(s) only found in [results UNNAMED] file: ", paste(extra_in_unamed, collapse = "\n\t\t - "))
        }
      }
    }
  }

  if(ic > 0) stop("\nPROBLEMS merging datasets")

  if(verbose) message("\n## MERGE")
  if(verbose) message("\nAll metabolomics datasets + basic phenotypic information")
  if(verbose) message("(it might take some time)")

  all_merged <- merge_all_metabolomics(m_m_n = m_m_n,
                                       m_m_u = m_m_u,
                                       m_s_n = m_s_n,
                                       r_n = r_m_n,
                                       r_u = r_m_u,
                                       phase = phase)

  return(all_merged)

}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Merge metabolomics metadata named and unnamed
#'
#' @description Merge metabolomics metadata
#' @param m_m_n (df) metabolomics metadata named
#' @param m_m_u (char) metabolomics metadata unnamed
#' @return (data.frame) merged metadata metabolites
#' @examples{
#' m_m <- merge_metabolomics_metadata(m_m_n = metadata_metabolites_named,
#'                                    m_m_u = metadata_metabolites_unnamed)
#' }
#' @export
merge_metabolomics_metadata <- function(m_m_n, m_m_u){

  colnames(m_m_n) <- tolower(colnames(m_m_n))

  m_m_n$id <- "named"
  m_m_u$id <- "unnamed"

  colnames(m_m_u) <- tolower(colnames(m_m_u))
  m_m_u$refmet_name <- NA
  m_m_u$formula <- NA
  m_m_u$formula <- NA
  m_m_u$neutral_mass <- NA

  right_columns <- c("metabolite_name", "refmet_name", "mz", "rt", "formula", "neutral_mass", "id")
  m_m_n <- m_m_n[right_columns]
  m_m_u <- m_m_u[right_columns]
  m_m <- rbind(m_m_n, m_m_u)

  m_m$metabolite <- ifelse(m_m$id == "named", m_m$refmet_name, m_m$metabolite_name)

  return(m_m)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Merge phenotypic and metabolics results
#'
#' @description Merge phenotypic data (phenotypes_pass1a06_short) and
#' metabolomics merged results and metadata
#' @param df_long (data.frame) Long format of a metabolomics merged results
#' @return (data.frame) Merged file, including the following columns:
#' \describe{
#'   \item{sample_id}{Sample Id, including vial_label and site specific QC ids}
#'   \item{sample_type}{Metabolomics sample types. Check metabolomics data transfer guidelines}
#'   \item{sample_order}{Order of injection on Mass Spec}
#'   \item{metabolite_name}{Given name by every lab}
#'   \item{refmet_name}{Map of the metabolite name to the Metabolomics RefMet database}
#'   \item{mz}{mass over charge}
#'   \item{rt}{retention time}
#'   \item{formula}{chemical formula}
#'   \item{neutral_mass}{neutral mass}
#'   \item{id}{type of metabolite identification: "named", "unnamed"}
#'   \item{metabolite}{Merge "refmet" for "named" metabolites and "metabolite_name" for "unnamed" metabolites}
#'   \item{quantification}{Untargeted: Peak area, Targeted: absolute concentration (check "experimentalDetails" for unit)}
#'   \item{tissue_code}{MoTrPAC tissue code}
#'   \item{tissue_name}{Tissue name}
#'   \item{group_time_point}{Intervention group (Exercise / Control) +  time point}
#'   \item{sex}{Animal Sex}
#'   \item{site_code}{Chemical Analysis Site (CAS) short abbreviation}
#'   \item{group}{Intervention group: Exercise / Control}
#'   \item{condition}{Sex + group + time-point}
#'   \item{bioreplicate}{Sex + group + time-point + sample_order}
#' }
#' @export
merge_phenotype_metabolomics <- function(df_long){

  df_merged <- merge(df_long, phenotypes_pass1a06_short, by.x = "sample_id", by.y = "vial_label", all.x = TRUE)
  cas <- unique(df_merged$site_code[!is.na(df_merged$site_code)])
  if(length(cas) != 1){
    stop("Problem with `cas_site` code from phenotypic data. Please, report this error to the BIC")
  }

  # Add labels, groups for plotting------
  df_merged$site_code <- cas
  df_merged$sex <- ifelse(is.na(df_merged$sex), "QC", df_merged$sex)
  df_merged$group <- ifelse(is.na(df_merged$sex), "QC", df_merged$group_time_point)
  df_merged$group <- ifelse(is.na(df_merged$group), "QC", df_merged$group)
  df_merged$group <- gsub("(Control)(.*)","\\1", df_merged$group)
  df_merged$group <- gsub("(Exercise)(.*)","\\1", df_merged$group)
  df_merged$group_time_point <- ifelse(is.na(df_merged$group_time_point), "QC", df_merged$group_time_point)
  df_merged$condition <- paste0(df_merged$sex, "_", df_merged$group_time_point)

  # Individual Controls
  df_merged$condition <- gsub("Female_Control - 00.0 hr", "F_Con_000h", df_merged$condition)
  df_merged$condition <- gsub("Female_Control - 07 hr", "F_Con_07h", df_merged$condition)
  # Combine controls
  # df_merged$condition <- gsub("Female_Control - 00.0 hr", "F_Con", df_merged$condition)
  # df_merged$condition <- gsub("Female_Control - 07 hr", "F_Con", df_merged$condition)

  df_merged$condition <- gsub("Female_Exercise - 00.0 hr", "F_Exe_000h", df_merged$condition)
  df_merged$condition <- gsub("Female_Exercise - 00.5 hr", "F_Exe_005h", df_merged$condition)
  df_merged$condition <- gsub("Female_Exercise - 01 hr", "F_Exe_01h", df_merged$condition)
  df_merged$condition <- gsub("Female_Exercise - 04 hr", "F_Exe_04h", df_merged$condition)
  df_merged$condition <- gsub("Female_Exercise - 07 hr", "F_Exe_07h", df_merged$condition)
  df_merged$condition <- gsub("Female_Exercise - 24 hr", "F_Exe_24h", df_merged$condition)
  df_merged$condition <- gsub("Female_Exercise - 48 hr", "F_Exe_48h", df_merged$condition)
  # Individual Controls
  df_merged$condition <- gsub("Male_Control - 00.0 hr", "M_Con_000h", df_merged$condition)
  df_merged$condition <- gsub("Male_Control - 07 hr", "M_Con_07h", df_merged$condition)
  # combined controls
  # df_merged$condition <- gsub("Male_Control - 00.0 hr", "M_Con", df_merged$condition)
  # df_merged$condition <- gsub("Male_Control - 07 hr", "M_Con", df_merged$condition)

  df_merged$condition <- gsub("Male_Exercise - 00.0 hr", "M_Exe_000h", df_merged$condition)
  df_merged$condition <- gsub("Male_Exercise - 00.5 hr", "M_Exe_005h", df_merged$condition)
  df_merged$condition <- gsub("Male_Exercise - 01 hr", "M_Exe_01h", df_merged$condition)
  df_merged$condition <- gsub("Male_Exercise - 04 hr", "M_Exe_04h", df_merged$condition)
  df_merged$condition <- gsub("Male_Exercise - 07 hr", "M_Exe_07h", df_merged$condition)
  df_merged$condition <- gsub("Male_Exercise - 24 hr", "M_Exe_24h", df_merged$condition)
  df_merged$condition <- gsub("Male_Exercise - 48 hr", "M_Exe_48h", df_merged$condition)
  df_merged$condition <- gsub("Male_Exercise - 48 hr", "M_Exe_48h", df_merged$condition)

  df_merged$condition <- gsub("QC_QC-Reference", "QC_Reference", df_merged$condition)
  df_merged$condition <- gsub("QC_QC-DriftCorrection", "QC_DriftCorrection", df_merged$condition)
  df_merged$condition <- gsub("QC_QC-Pooled", "QC_Pooled", df_merged$condition)

  df_merged$bioreplicate <- paste0(df_merged$condition, "-", df_merged$sample_order)

  df_merged <- df_merged[,c("site_code",
                            "sample_id", "sample_type", "sample_order",
                            "tissue_code", "tissue_name",
                            "group_time_point", "sex",
                            "group", "condition", "bioreplicate",
                            "metabolite_name", "refmet_name", "mz", "rt", "formula", "neutral_mass",
                            "id", "quantification")]

  return(df_merged)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Merge all metabolomics files
#'
#' @description Merge all metabolomics datasets, including "results" and
#' "metadata" files, for both targeted and untargeted datasets
#' @param m_m_n (metabolomics metadata named)
#' @param m_m_u (metabolomics metadata unnamed)
#' @param m_s_n (metabolomics sample named)
#' @param r_n (results named)
#' @param r_u (results unnamed)
#' @param phase (MoTrPAC Animal phase. Eg. PASS1A-06)
#' @return (data.frame) Merged data frame long format
#' @examples
#' plasma.untargeted.merged <- merge_all_metabolomics(
#'        m_m_n = metadata_metabolites_named,
#'        m_m_u = metadata_metabolites_unnamed,
#'        m_s_n = metadata_sample_named,
#'        r_n = results_named,
#'        r_u = results_unnamed,
#'        phase = "PASS1A-06")
#' @export
merge_all_metabolomics <- function(m_m_n,
                                   m_m_u = NULL,
                                   m_s_n,
                                   r_n,
                                   r_u = NULL,
                                   phase){

  # # Debug
  # m_m_n = metadata_metabolites_named
  # m_m_u = metadata_metabolites_unnamed
  # m_s_n = metadata_sample_named
  # r_n = results_named
  # r_u = results_unnamed

  raw_file = NULL

  # metabolites_metadata----
  if(is.null(m_m_u)){
    right_columns <- c("metabolite_name", "refmet_name", "mz", "rt", "formula", "neutral_mass", "id")
    m_m_n$id <- "named"
    m_m <- m_m_n[right_columns]
  }else{
    m_m <- merge_metabolomics_metadata(m_m_n, m_m_u)
  }

  # merge results metabolites
  if(is.null(r_u)){
    r_m <- r_n
  }else{
    r_m <- rbind(r_n, r_u[names(r_n)])
  }

  r_long <- tidyr::pivot_longer(r_m,
                                cols = -c("metabolite_name"),
                                names_to = "sample_id",
                                values_to = "quantification",
                                values_drop_na = FALSE)

  # Merge with metadata metabolites-----
  before <- dim(r_long)[1]
  r_long <- merge(m_m, r_long, by = "metabolite_name")
  after <- dim(r_long)[1]
  if(before != after) stop("PROBLEMS WITH DIMENSIONS after merging!")

  # Merge with sample metadata-----
  # Remove raw file
  m_s_n <- subset(m_s_n, select = -raw_file)

  before <- dim(r_long)[1]
  r_long <- merge(m_s_n, r_long, by = "sample_id")
  after <- dim(r_long)[1]
  if(before != after) stop("PROBLEMS WITH DIMENSIONS after merging!")

  if(phase == "PASS1A-06"){
    r_long_pheno <- merge_phenotype_metabolomics(df_long = r_long)
  }else{
    r_long_pheno <- r_long
    message("\n(-) Phenotypic data not available yet in this package for ", phase)
  }

  return(r_long_pheno)
}





