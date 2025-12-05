
#' @title Check olink metadata proteins file
#'
#' @description check whether metadata-proteins.txt file is following guidelines
#' @param df (data.frame) metadata_proteins df
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param validate_uniprot (logical) if `TRUE`, check if all uniprot ids are valid
#' connecting to Uniprot. Note: depending on the number of ids, it might take
#' several finish to complete the validation
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' \dontrun{
#' check_metadata_proteins(df = metadata_proteins)
#' }
#' }
#' @export
check_metadata_proteins <- function(df,
                                    return_n_issues = FALSE,
                                    validate_uniprot = FALSE,
                                    verbose = TRUE){
  
  # issue_count
  ic <- 0
  
  df <- filter_required_columns(df = df,
                                type = "olproteins",
                                verbose = verbose)
  
  
  # Evaluate every column
  
  if( "olink_id" %in% colnames(df) ){
    if( length(unique(df$olink_id)) != dim(df)[1] ){
      duplis_details <- df$olink_id[duplicated(df$olink_id)]
      duplis <- length(unique(duplis_details))
      if(verbose) message("   - (-) `olink_id` non-unique values detected: ", duplis)
      if(verbose) message("\t\t - ", paste(unique(duplis_details), collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("  + (+) `olink_id`: unique values: OK")
    }
    
    if( check_missing_values(df, "olink_id") ){
      if(verbose) message("   - (-) `olink_id`: NA values detected: FAIL")
      ic <- ic + 1
    }
  }else{
    if(verbose) message("   - (-) `olink_id`: is missed: FAIL")
    ic <- ic + 1
  }
  
  # uniprot_entry
  if("uniprot_entry" %in% colnames(df)){
    if( length(unique(df$uniprot_entry)) != dim(df)[1] ){
      duplis_details <- df$uniprot_entry[duplicated(df$uniprot_entry)]
      duplis <- length(unique(duplis_details))
      if(verbose) message("   - ( ) `uniprot_entry` non-unique values detected (n duplications = ", duplis, "). This is OK")
      if(verbose) message("\t\t - ", paste(unique(duplis_details), collapse = "\n\t\t - "))
    }else{
      if(verbose) message("  + (+) `uniprot_entry` unique values: OK")
    }
    
    if(validate_uniprot){
      if(verbose) message("  + Validating every `uniprot_entry`. Connecting to Uniprot database (please, be patient as this step might take several minutes)")
      all_true <- validate_uniprot_ids_with_uniprot(ids = unique(df$uniprot_entry))
      if(all_true){
        if(verbose) message("  + (+) `uniprot_entry` ids found in Uniprot database: OK")
      } else{
        if(verbose) message("  + (-) One or many `uniprot_entry` ids not found in Uniprot database: FAIL")
        ic <- ic + 1
      }
    }
  }else{
    if(verbose) message("   - (-) `uniprot_entry` column missed: FAIL")
    ic <- ic + 1
  }
  
  # assay
  if("assay" %in% colnames(df)){
    
    if(length(unique(df$assay)) != dim(df)[1]){
      
      duplis_details <- df$assay[duplicated(df$assay)]
      duplis <- length(unique(duplis_details))
      
      if(verbose) message("   - ( ) `assay` non-unique values detected (n duplications = ", duplis, "). This is OK")
      if(verbose) message("\t - ", paste(unique(duplis_details), collapse = "\n\t - "))
      
    }else{
      if(verbose) message("  + (+) `assay` unique values: OK")
    }
    
    if( check_missing_values(df, "assay") ){
      if(verbose) message("   - (-) `assay` NA values detected: FAIL")
      ic <- ic + 1
    }
  }else{
    if(verbose) message("   - (-) `assay` column missed: FAIL")
    ic <- ic + 1
  }
  
  if("missing_freq" %in% colnames(df)){
    if( !all(is.numeric(df$missing_freq)) ){
      if(verbose) message("   - (-) `missing_freq` non numeric values found: FAIL")
      nonnum <- df$missing_freq[which(!grepl('^[0-9]', df$missing_freq))]
      if(verbose) message("\n\t - ", paste(nonnum, collapse = "\n\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("  + (+) `missing_freq` all numeric: OK")
    }
    
    if( check_missing_values(df, "missing_freq") ){
      if(verbose) message("   - (-) `missing_freq` NA values detected: FAIL")
      ic <- ic + 1
    }
    
  }else{
    if(verbose) message("   - (-) `missing_freq` column missed: FAIL")
    ic <- ic + 1
  }
  
  if("panel_name" %in% colnames(df)){
    if(verbose){
      message("  + (+) `panel_name` checking available panels:")
      message("\t - ", paste(unique(df$panel_name), collapse = "\n\t - "))
    } 
    
    if( check_missing_values(df, "panel_name") ){
      if(verbose) message("   - (-) `panel_name` NA values detected: FAIL")
      ic <- ic + 1
    }
  }else{
    if(verbose) message("   - (-) `panel_name` column missed: FAIL")
    ic <- ic + 1
  }

  if("panel_lot_nr" %in% colnames(df)){
    if(verbose){
      message("  + (+) `panel_lot_nr` checking available panels:")
      message("\t - ", paste(unique(df$panel_lot_nr), collapse = "\n\t - "))
    } 
    
    if( check_missing_values(df, "panel_lot_nr") ){
      if(verbose) message("   - (-) `panel_lot_nr` NA values detected: FAIL")
      ic <- ic + 1
    }
  }else{
    if(verbose) message("   - (-) `panel_lot_nr` column missed: FAIL")
    ic <- ic + 1
  }
  
  if("normalization" %in% colnames(df)){
    if(verbose){
      message("  + (+) `normalization` checking available panels:")
      message("\t - ", paste(unique(df$normalization), collapse = "\n\t - "))
    } 
    
    if( check_missing_values(df, "normalization") ){
      if(verbose) message("   - (-) `normalization` NA values detected: FAIL")
      ic <- ic + 1
    }
  }else{
    if(verbose) message("   - (-) `normalization` column missed: FAIL")
    ic <- ic + 1
  }
  
  if(return_n_issues) return(ic)
  
} #end check_metadata_proteins


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title check olink metadata samples file
#'
#' @description check whether metadata_sample is following guidelines
#' @param df (data.frame) olink metadata samples ata
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' \dontrun{
#' check_metadata_samples(df = metadata_sample_named)
#' }
#' }
#' @export
check_metadata_samples_olink <- function(df,
                                         return_n_issues = FALSE,
                                         verbose = TRUE){
  
  olink_id = sample_order = unique_sample_order = plate_id = NULL
  
  # issue_count
  ic <- 0
  
  # filter only expected columns
  df <- filter_required_columns(df = df,
                                type = "olsamples",
                                verbose = verbose)
  
  if( "sample_id" %in% colnames(df) ){
    if( length(unique(df$sample_id)) != dim(df)[1] ){
      if(verbose) message("   - (-) `sample_id`: Non-unique values detected: FAIL")
      ic <- ic + 1
    }else{
      if(verbose) message("  + (+) `sample_id` seems OK")
    }
  }else{
    if(verbose) message("   - (-) `sample_id` is missed: FAIL")
    ic <- ic + 1
  }
  
  # sample_type: st
  esample_types <- c("Sample", "QC-Pooled", "QC-Reference", "QC-Blank",
                     "QC-Identification", "QC-InternalStandard", "QC-PreRun",
                     "QC-ExternalStandard", "QC-DriftCorrection", "QC-ReCAS", 
                     "QC-PlateControl")
  
  if("sample_type" %in% colnames(df)){
    if(!all(df$sample_type %in% esample_types)){
      if(verbose) message("   - (-) Error: undefined sample types: ", appendLF = FALSE)
      if(verbose) message("\n\t\t - ", paste(setdiff(df$sample_type, esample_types), collapse = "\n\t\t - "))
      ic <- ic + 1
    }else{
      if(verbose) message("  + (+) `sample_type` seems OK")
    }
  }else{
    if(verbose) message("   - (-) `sample_type` column missed: FAIL")
    ic <- ic + 1
  }
  
  if( "plate_id" %in% colnames(df) ){
    if(verbose) message("  + (+) `plate_id` is available: OK")
    if( "sample_order" %in% colnames(df) ){
      if(!all(is.numeric(df$sample_order))){
        if(verbose) message("   - (-) `sample_order` non numeric values found: ", appendLF = FALSE)
        nonnum <- df$sample_order[which(!grepl('^[0-9]', df$sample_order))]
        if(verbose) message("\t\t - ", paste(nonnum, collapse = "\n\t\t - "))
        ic <- ic + 1
      }else{
        if(verbose) message("  + (+) `sample_order` is numeric: OK")
        
        # Check if the sample order is unique for each plate-id
        non_unique_results <- df %>%
          group_by(plate_id) %>%
          summarise(unique_sample_order = all(length(unique(sample_order)) == length(sample_order))) %>%
          filter(unique_sample_order == FALSE)
        
        # Print plate_ids with non-unique sample_order values
        if (nrow(non_unique_results) > 0) {
          if(verbose) message("   - (-) `plate_id` with non-unique `sample_order` values: FAIL")
          ic <- ic + 1
          if(verbose) message("\t\t - ", paste(non_unique_results$plate_id, collapse = "\n\t\t - "))
        } else {
          if(verbose) message("  + (+) All `plate_id` values have unique `sample_order` values: OK\n")
        }
      }
    }else{
      if(verbose) message("   - (-) `sample_order` column missed: FAIL")
      ic <- ic + 1
    }
    

  }else{
    if(verbose) message("   - (-) `plate_id` column missed: FAIL")
    ic <- ic + 1
  }
  
  if(return_n_issues) return(ic)
} #end check_metadata_sample_olink


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Check results file for assays
#'
#' @description Checks whether the results file follows guidelines for Olink and IMM assays.
#' @param df (data.frame) The results data frame to check.
#' @param assay_type (character) The type of assay, either `"olink"` or `"imm"`.
#' @param return_n_issues (logical) If `TRUE`, returns the number of issues identified.
#' @param verbose (logical) If `TRUE` (default), displays messages during the checking process.
#' @return (int) Number of issues identified if `return_n_issues` is `TRUE`.
#' @examples
#' \dontrun{
#' check_results_assays(df = results_df, assay_type = "olink")
#' check_results_assays(df = results_df, assay_type = "imm")
#' }
#' @export
check_results_assays <- function(df,
                          assay_type = c("olink", "lab"),
                          return_n_issues = FALSE,
                          verbose = TRUE) {
  
  # Match the assay_type argument
  assay_type <- match.arg(assay_type)
  
  # Initialize issue count
  ic <- 0
  
  is_empty_df <- ncol(df) == 0 && nrow(df) == 0
  
  if (is_empty_df) {
    if (verbose) message("   - (-) `results` data frame is empty: FAIL")
    if (return_n_issues) {
      ic <- 10
      return(ic)
    } else {
      return("Data frame is empty")
    }
  }
  
  # Determine the identifier column based on assay type
  id_column <- if (assay_type == "olink") {
    "olink_id"
  } else if (assay_type == "lab") {
    "analyte_name"
  }
  
  # Check for identifier column
  if (id_column %in% colnames(df)) {
    if (length(unique(df[[id_column]])) != nrow(df)) {
      if (verbose) message("   - (-) `", id_column, "`: Non-unique values detected: FAIL")
      ic <- ic + 1
    } else {
      if (verbose) message("  + (+) `", id_column, "` unique values: OK")
    }
  } else {
    if (verbose) message("   - (-) `", id_column, "` column missing: FAIL")
    ic <- ic + 1
  }
  
  # Get the names of all columns except the identifier column
  columns_to_check <- setdiff(names(df), id_column)
  
  # Check if each of these columns is numeric
  are_numeric <- sapply(df[columns_to_check], is.numeric)
  
  # Print results
  if (all(are_numeric)) {
    if (verbose) message("  + (+) All measurement columns are numeric: OK")
  } else {
    if (verbose) message("   - (-) `results` contains non-numeric columns: FAIL")
    non_numeric_columns <- names(df[columns_to_check])[!are_numeric]
    if (verbose) message("\t\t - ", paste(non_numeric_columns, collapse = "\n\t\t - "))
    ic <- ic + 1
  }
  
  # Calculate counts
  na_counts <- sapply(df[columns_to_check], function(column) sum(is.na(column)))
  zero_counts <- sapply(df[columns_to_check], function(column) sum(column == 0, na.rm = TRUE))
  total_na_count <- sum(na_counts)
  total_zero_count <- sum(zero_counts)
  total_values_count <- sum(sapply(df[columns_to_check], function(column) length(column))) - total_na_count
  
  # Print the result
  if (verbose) {
    message(paste0("  + ( ) Number of zeros in dataset: ", total_zero_count, " (out of ", total_values_count, " values)"))
    message(paste0("  + ( ) Number of NAs in dataset: ", total_na_count, " (out of ", total_values_count + total_na_count, " values)"))
  }
  
  if (return_n_issues) return(ic)
} # End of check_results_assays


#' @title Load and Process Olink Batch Data
#'
#' @description 
#' This function loads Olink batch data from the specified input directory.
#' It performs quality checks on the data and loads specific files related to
#' Olink data, including metadata for proteins and samples, and the results file.
#' It also integrates validation checks and warns if there are too many issues
#' identified in the data.
#'
#' @param input_results_folder A string representing the path to the folder
#' containing Olink batch data to be loaded and processed.
#' @param verbose Logical; if `TRUE`, prints detailed messages 
#' during the loading process.
#'
#' @return A list containing data frames for metadata of proteins (m_p),
#' metadata of samples (m_s), and Olink results (r_o). If certain files
#' are not available, the corresponding entries in the list will be NULL.
#' @examples
#' \dontrun{
#' list_of_df <- load_olink_batch(input_results_folder = "/path/to/PROCESSED_YYYYMMDD/")
#' }
#' @export
load_olink_batch <- function(input_results_folder,
                             verbose = TRUE){
                                    
  
  m_p = m_s = r_o = NULL
  
  # Validations----
  phase <- validate_phase(input_results_folder)
  processfolder <- validate_processFolder(input_results_folder)
  assay <- validate_assay(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)
  
  total_issues <- validate_olink(input_results_folder = input_results_folder, 
                                 cas = "broad_rg", 
                                 return_n_issues = TRUE, 
                                 verbose = verbose)
  
  if(total_issues > 0){
    message("\n\tWARNING!!! Too many issues identified (", total_issues,"). 
            This batch should not be processed until the issues are solved")
  }
  
  vial_label <- NA
  qc_samples <- NA
  
  # Load olink----
  if(verbose) message("# LOAD OLINK BATCH")
  if(verbose) message("+ Site: Broad, Gerszten Lab")
  if(verbose) message("+ Folder: `", paste0(input_results_folder),"`")
  
# qc metadata-proteins----
  if(verbose) message("\n## QC `metadata_proteins`\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata-proteins.txt",
                     verbose = verbose)
  f_mp <- lista$flag
  if(f_mp){
    m_p_f <-lista$filename
    m_p <- lista$df
                                           
  }else{
    if(verbose) message("   - (-) `metadata_proteins.txt` file not available")
    ic_m_p <- 10
  }
  
  # qc metadata-samples------
  if(verbose) message("\n## QC `metadata-samples.txt`\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata-samples.txt",
                     verbose = verbose)
  f_ms <- lista$flag
  if(f_ms){
    m_s_f <-lista$filename
    m_s <- lista$df
  }else{
    if(verbose) message("   - (-) `metadata-samples.txt` file not available")
    ic_m_s <- 10
  }
  
  # qc results olink------
  if(verbose) message("\n## QC `results.txt`\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results.txt",
                     verbose = verbose)
  f_r <- lista$flag
  if( f_r ){
    r_o_f <-lista$filename
    r_o <- lista$df
  }else{
    if(verbose) message("   - (-) `results.txt` file not available")
    ic_r <- 10
  }
  
  # RETURN list of dfs-----
  
  listdf <- list ("m_p" = m_p,
                  "m_s" = m_s,
                  "r_o" = r_o)
  
  return(listdf)
}


#' Validate Olink Data
#'
#' This function validates Olink data based on a set of criteria, including
#' folder structure, metadata, and file contents. It supports optional
#' functionalities like creating a PDF report and DMAQC validation 
#' (data only available at the BIC)
#'
#' @param input_results_folder A string representing the path to the folder
#'        containing Olink results to be validated.
#' @param cas A character string indicating the CAS number.
#' @param return_n_issues Logical; if `TRUE`, the function returns the number of
#'        detected issues.
#' @param full_report Logical; if `TRUE`, generates a full report of the 
#'        validation process.
#' @param f_proof Logical; if `TRUE`, generates proof plots for data validation.
#' @param printPDF Logical; if `TRUE` and `f_proof` is `TRUE`, saves the plots
#' to a PDF file, and in such case, then provide the desired path to output
#' the PDF file in the argument `out_qc_folder`
#' @param out_qc_folder Optional; a string specifying the path to the folder
#' where output PDF should be saved (only if `printPDF = TRUE`). 
#' Default: current working directory
#' @param dmaqc_shipping_info (char) File path to the DMAQC file. 
#' Only the BIC can use this argument
#' @param dmaqc_phase2validate (char) Provide phase to validate. This argument
#' is not required since it should be extracted from the input folder or from the 
#' new required file `metadata_phase.txt`. Please, ignore. 
#' However, if this argument is provided,
#' it will take priority and this will be the phase.
#' @param validate_uniprot Logical; if `TRUE`, validates against the UniProt
#'        database.
#' @param verbose Logical; if `TRUE`, prints detailed messages during validation.
#'
#' @return Depending on the settings, this function may return the number of
#' issues found, generate reports or plots, or simply perform the 
#' validation without returning anything.
#' @examples
#' \dontrun{
#' validate_olink("/path/to/results", cas = "broad_rg", return_n_issues = TRUE)
#' }
#' @export
validate_olink <- function(input_results_folder,
                           cas,
                           return_n_issues = FALSE,
                           full_report = FALSE,
                           f_proof = FALSE,
                           printPDF = FALSE,
                           out_qc_folder = NULL,
                           dmaqc_shipping_info = NULL,
                           dmaqc_phase2validate = FALSE,
                           validate_uniprot = FALSE,
                           verbose = TRUE){
                                  
  
  olink_id = sample_order = plate_id = sample_id = NULL
  
  # validate folder structure -----
  validate_cas(cas = cas)
  processfolder <- validate_processFolder(input_results_folder)
  assay <- validate_assay(input_results_folder)
  phase <- validate_phase(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)
  batch_folder <- validate_batch(input_results_folder)
  
  # issue_count-----
  ic <- 0
  ic_m_p <- 0 # ic for metadata protein file
  ic_m_s <- 0 # ic for metadata sample file
  ic_r <- 0
  ic_man <- 0 # required manifest f
  
  if(is.null(dmaqc_shipping_info)){
    ic_vl <- "missed"
  } else{
    ic_vl <- NA
  }
  
  vial_label <- NA
  qc_samples <- NA
  
  input_results_folder <- normalizePath(input_results_folder)
  
  input_folder_short <- regmatches(input_results_folder, regexpr("(HUMAN|PASS).*RESULTS_[0-9]{8}", input_results_folder))
  if(purrr::is_empty(input_folder_short)){
    if(verbose) message("\nThe PROCESSED_YYYYMMDD folder full path is not correct. Example:")
    if(verbose) message("/full/path/to/folder/HUMAN-PRECOVID/T02/BATCH1_20190822/RESULTS_20200302")
    stop("Input folder not according to guidelines")
  }
  
  if(verbose) message("# OLINK QC report\n\n")
  if(verbose) message("+ Site: ", cas)
  if(verbose) message("+ Folder: `",paste0(input_folder_short),"`")
  
  is_mp <- check_metadata_phase_file(input_results_folder = input_results_folder, 
                                     verbose = verbose)
  if(!is_mp){
    ic <- ic + 1
  }
  
  # Set phase-----
  dmaqc_phase2validate <- set_phase(input_results_folder = input_results_folder, 
                                    dmaqc_phase2validate = dmaqc_phase2validate,
                                    verbose = verbose)
  
  phase2file <- generate_phase_details(phase_metadata = dmaqc_phase2validate)
  
  # qc metadata-proteins----
  if(verbose) message("\n## QC `metadata_proteins`\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata-proteins.txt",
                     verbose = verbose)
  f_mp <- lista$flag
  if(f_mp){
    m_p_f <-lista$filename
    m_p <- lista$df
    ic_m_p <- check_metadata_proteins(df = m_p,
                                      return_n_issues = TRUE,
                                      validate_uniprot = validate_uniprot,
                                      verbose = verbose)
                                           
  }else{
    if(verbose) message("   - (-) `metadata_proteins.txt` file not available")
    ic_m_p <- 10
  }
  
  # qc metadata-samples------
  if(verbose) message("\n## QC `metadata-samples.txt`\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "metadata-samples.txt",
                     verbose = verbose)
  f_ms <- lista$flag
  if(f_ms){
    m_s_f <-lista$filename
    m_s <- lista$df
    ic_m_s <- check_metadata_samples_olink(df = m_s,
                                           return_n_issues = TRUE,
                                           verbose = verbose)
    # Extract the number of samples
    if(!is.null(m_s)){
      #Double check that the columns are there
      if( all(c("sample_id", "sample_type") %in% colnames(m_s)) ){
        vial_label <- length(m_s$sample_id[which(m_s$sample_type == "Sample")])
        qc_samples <- length(m_s$sample_id[which(m_s$sample_type != "Sample")])
      }
    }
  }else{
    if(verbose) message("   - (-) `metadata-samples.txt` file not available")
    ic_m_s <- 10
  }
  
  # qc results olink----
  if(verbose) message("\n## QC `results.txt`\n")
  
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results.txt",
                     verbose = verbose)
  f_r <- lista$flag
  if( f_r ){
    r_o_f <-lista$filename
    r_o <- lista$df
    ic_r <- check_results_assays(df = r_o,
                                return_n_issues = TRUE,
                                assay_type = "olink",
                                verbose = verbose)
    
  }else{
    if(verbose) message("   - (-) `results.txt` file not available")
    ic_r <- 10
  }
  
  # qc cross validation-----
  if(verbose) message("\n## Cross File Validation\n")
  if(ic_m_p == 0 && ic_m_s == 0 && ic_r == 0) {
    ic_c_f_v <- check_crossfile_validation(r_o = r_o, 
                                           m_s = m_s, 
                                           m_p = m_p, 
                                           assay_type = "olink",
                                           return_n_issues = TRUE, 
                                           verbose = verbose)
  }else{
    if(verbose) message("   - (-) File cross validation is not possible. PLease, fix issues affecting any of the files and run the validation again.")
  }
  
  # QC PLOTS------
  
  if(f_proof){
    
    if(verbose) message("\n\n## QC Plots\n")
    
    output_prefix <- paste0(cas, ".", tolower(phase2file), ".", tissue_code, ".",tolower(assay), ".", tolower(processfolder))
    output_prefix <- gsub("\\|", "_", output_prefix)
    
    if(f_r & f_ms & f_mp){
      # Ensure there are enough compounds to generate plots
      if( dim(r_o)[1] > 1 ){
        results_long <- r_o %>% tidyr::pivot_longer(cols = -c(olink_id),
                                                        names_to = "sample_id",
                                                        values_to = "value")
        
        results_long <- merge(m_s, results_long, by = c("sample_id"))
        results_long$sample_id <- as.character(results_long$sample_id)
        results_long$sample_id <- as.factor(results_long$sample_id)
        results_long$sample_type <- as.factor(results_long$sample_type)
        results_long <- results_long[which(results_long$value != 0),]
        results_long <- results_long[!is.na(results_long$value),]
        
        # sort by plate_id and sample_order
        results_long <- results_long %>%
          arrange(plate_id, sample_order) %>%
          mutate(sample_id_ordered = factor(sample_id, levels = unique(sample_id)))
        
        r_p <- merge(m_p, r_o, by = "olink_id")
        
        plot_basic_olink_qc(results = r_p, 
                            results_long = results_long,
                            out_qc_folder = out_qc_folder, 
                            output_prefix = output_prefix,
                            printPDF = printPDF,
                            verbose = verbose)
                                   
        
      }else{
        message("  (-) QC plots are not possible: not enough compounds")
      }
      
    }else{
      if(verbose) message("\n- (-) QC plots are not possible: critical datasets are missed")
    }
  } # qc plots
  
  # MANIFEST file-----
  
  if(verbose) message("\n## QC `file_manifest_YYYYMMDD.csv` (required)\n")
  
  batch <- gsub("(.*)(RESULTS.*)", "\\1", input_results_folder)
  
  file_manifest <- list.files(normalizePath(batch),
                              pattern="file_manifest.*csv",
                              ignore.case = TRUE,
                              full.names=TRUE,
                              recursive = TRUE)
  
  
  if(length(file_manifest) == 0){
    f_man <- FALSE
    if(verbose) message("   - (-) `file_manifest_YYYYMMDD.csv` file not available")
    ic <- ic + 1
  }else if(length(file_manifest) >= 1){
    first_manifest <- sort(basename(file_manifest), decreasing = TRUE)[1]
    file_manifest <- file_manifest[grep(first_manifest, file_manifest)]
    f_man <- TRUE
  }
  
  if(f_man){
    manifest <- read.csv(file_manifest)
    mani_columns <- c("file_name", "md5")
    if( all( mani_columns %in% colnames(manifest)) ){
      if(verbose) message("  + (+) `file_name, md5` columns available in manifest file")
      
      # Replace windows based backlash
      if(any(grepl("\\\\", manifest$file_name))){
        manifest$file_name <- gsub("\\\\", "/", manifest$file_name)
      }
      manifest$file_base <- basename(manifest$file_name)
      
      if(f_mp){
        metadata_proteins_file <- basename(manifest$file_name[grepl("metadata-proteins", manifest$file_name)])[1]
        tocheck <- basename(m_p_f)
        if(!is.na(metadata_proteins_file)){
          if( tocheck == metadata_proteins_file){
            if(verbose) message("  + (+) `metadata-proteins` file included in manifest: OK")
          }else{
            if(verbose) message("   - (-) `metadata-proteins` file  is not included in manifest file: FAIL")
            ic_man <- ic_man + 1
          }
        }else{
          if(verbose) message("   - (-) `metadata-proteins` file  is not included in manifest file: FAIL")
          ic_man <- ic_man + 1
        }
      }
      
      if(f_ms){
        metadata_samples_file <- basename(manifest$file_name[grepl("metadata-samples", manifest$file_name)])[1]
        tocheck <- basename(m_s_f)
        if(!is.na(metadata_samples_file)){
          if( tocheck == metadata_samples_file ){
            if(verbose) message("  + (+) `metadata-samples` file included in manifest: OK")
          }else{
            if(verbose) message("   - (-) `metadata-samples` file is not included in manifest file: FAIL")
            ic_man <- ic_man + 1
          }
        }else{
          if(verbose) message("   - (-) `metadata-samples` file is not included in manifest file: FAIL")
          ic_man <- ic_man + 1
        }
      } # f_msn
      
      if(f_r){
        results_file <- basename(manifest$file_name[grepl("results", manifest$file_name)])[1]
        tocheck <- basename(r_o_f)
        if(!is.na(results_file)){
          if( tocheck == results_file ){
            if(verbose) message("  + (+) `results` file included in manifest: OK")
          }else{
            if(verbose) message("   - (-) `results` file is not included in manifest file: FAIL")
            ic_man <- ic_man + 1
          }
        }
        
      }
      
      if( any(is.na(manifest$md5)) ){
        if(verbose) message("   - (-) MD5 column contains NA values: FAIL")
        ic_man <- ic_man + 1
      }
      
    }else{
      if(verbose) message("   - (-) Not all the required columns are available: FAIL")
      ic_man <- ic_man + 1
    }
  }else{
    if(verbose) message("   - (-) MANIFEST (REQUIRED) FILE NOT FOUND (`file_manifest_DATE.txt`). Please, check guidelines")
    ic_man <- ic_man + 6
    ic <- ic + 1
  } # QC manifest
  
  # DMAQC validation -----
  if(verbose) message("\n\n## DMAQC validation\n")
  failed_samples <- check_failedsamples(input_results_folder = input_results_folder, 
                                        verbose = verbose)
  
  # Validate vial labels from DMAQC
  if( is.na(ic_vl) ){
    if(f_ms){
      vl_results <- m_s$sample_id[which(m_s$sample_type == "Sample")]
      outfile_missed_viallabels <- paste0(cas, ".", tolower(phase2file), ".", tissue_code, ".",tolower(assay), ".", tolower(processfolder))
      outfile_missed_viallabels <- gsub("\\|", "_", outfile_missed_viallabels)
      ic_vl <- check_viallabel_dmaqc(vl_submitted = vl_results,
                                     tissue_code = tissue_code,
                                     cas = cas,
                                     phase = dmaqc_phase2validate,
                                     failed_samples = failed_samples,
                                     dmaqc_shipping_info = dmaqc_shipping_info,
                                     out_qc_folder = out_qc_folder,
                                     outfile_missed_viallabels = outfile_missed_viallabels,
                                     return_n_issues = TRUE,
                                     verbose = verbose)
    }
  }
  
  
  # RETURN report-----
  if(ic > 4){
    message("\nTOTAL NUMBER OF CRITICAL ERROR: ", ic,"\n")
    message("WARNING: Too many errors. Revise input folder")
  }
  
  batchversion <- stringr::str_extract(string = input_results_folder, pattern = "BATCH.*_[0-9]+/RESULTS_[0-9]+")
  
  qc_date <- format(Sys.time(), "%Y%m%d")
  t_name <- bic_animal_tissue_code$bic_tissue_name[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)]
  
  if(return_n_issues){
    total_issues <- sum(ic, ic_man, ic_m_p, ic_m_s, ic_r, na.rm = TRUE)

    if(verbose) message("\nTOTAL NUMBER OF ISSUES: ", total_issues,"\n")
    if(full_report){
      reports <- data.frame(cas = cas,
                            phase = dmaqc_phase2validate,
                            tissue = tissue_code,
                            t_name = t_name,
                            assay = assay,
                            version = batchversion,
                            vial_label = vial_label,
                            qc_samples = qc_samples,
                            dmaqc_valid = ic_vl,
                            critical_issues = ic,
                            manifest = ic_man,
                            m_prot = ic_m_p,
                            m_sample = ic_m_s,
                            results = ic_r,
                            qc_date = qc_date)
      return(reports)
    }else{
      return(total_issues)
    }
  }
} # Validate OLINK


#' @title Write olink data release
#'
#' @description Write out olink data releases. Doesn't check whether
#' data has been submited according to guidelines
#' @param input_results_folder (char) Path to the RESULTS_YYYYMMDD folder
#' @param folder_name (char) output folder name.
#' @param folder_root (char) absolute path to write the output folder. 
#' Default: current directory
#' @param version_file (char) file version number (v#.#)
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return bic release folder/file structure 
#' `PHASE/OMICS/TCODE_NAME/ASSAY/` and file names, including:
#' - `motrpac_phase-code_tissuecode_assay_file-details-version.txt`
#'  where files-details can be:
#' - `metadata-proteins`, 
#' - `metadata-samples`,
#' - `results`
#' @examples
#' \dontrun{
#' write_olink_releases(
#'    input_results_folder = "/full/path/to/RESULTS_YYYYMMDD/")
#' }
#' @export
write_olink_releases <- function(input_results_folder,
                                 folder_name = "motrpac_release",
                                 folder_root = NULL,
                                 version_file = "v1.0",
                                 verbose = TRUE){
  
  # Get names from input_results_folder------
  assay <- validate_assay(input_results_folder)
  
  phase <- validate_phase(input_results_folder)
  
  phase_metadata <- set_phase(input_results_folder = input_results_folder, 
                              dmaqc_phase2validate = FALSE, 
                              verbose = FALSE)
  phase_details <- generate_phase_details(phase_metadata)
  
  tissue_code <- validate_tissue(input_results_folder)
  
  folder_tissue <- bic_animal_tissue_code$tissue_name_release[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)]
  
  # # or make a function:
  # get_folder_tissue <- function(tissue_code) {
  #   folder_tissue <- MotrpacBicQC::bic_animal_tissue_code$tissue_name_release[
  #     which(MotrpacBicQC::bic_animal_tissue_code$bic_tissue_code == tissue_code)
  #   ]
  #   return(folder_tissue)
  # }
  # 
  # folder_tissue2 <- get_folder_tissue(tissue_code)
  
  if(length(assay_codes$assay_code[which(assay_codes$submission_code == assay)]) == 1 ){
    folder_assay <- assay_codes$assay_code[which(assay_codes$submission_code == assay)]
  }else{
    stop("ASSAY code ", assay, " not available in `assay_codes`")
  }
  
  if(verbose) message("+ Writing out: ", phase_details, " ", tissue_code, " ", assay, " files", appendLF = FALSE)
  
  # Load olink datasets----
  olink_df <- load_olink_batch(input_results_folder = input_results_folder,
                                verbose = FALSE)
  
  # Create output folder-------
  if (is.null(folder_root)){
    folder_root <- getwd()
  }else{
    folder_root <- create_folder(folder_root)
  }
  
  # Exception for PASS1C-06: the main folder is pass1a
  if(phase_details == "pass1c-06"){
    phase_folder_release <- "pass1a-06"
  }else{
    phase_folder_release <- phase_details
  }
  
  output_folder <- file.path(folder_root, folder_name, phase_folder_release, "proteomics-targeted", folder_tissue, folder_assay)
  
  if(!dir.exists(file.path(output_folder))){
    dir.create(file.path(output_folder), recursive = TRUE)
  }
  
  file_name_shared <- paste0("motrpac_",
                             phase_details, "_",
                             folder_tissue, "_",
                             folder_assay)
  
  # Create and write FILES-----
  metadata_proteins <- file.path(output_folder, paste0(file_name_shared,"_metadata-proteins_", version_file, ".txt"))
  metadata_samples <- file.path(output_folder, paste0(file_name_shared,"_metadata-samples_", version_file, ".txt"))
  results <- file.path(output_folder, paste0(file_name_shared,"_results_", version_file, ".txt"))
  
  write.table(olink_df$m_p, metadata_proteins, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(olink_df$m_s, metadata_samples, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(olink_df$r_o, results, row.names = FALSE, sep = "\t", quote = FALSE)
  
  if(verbose) message("...done!")
}



