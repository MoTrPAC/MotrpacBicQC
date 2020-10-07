
utils::globalVariables(
  c("col_name",
    "confident_score",
    "confident_site",
    "entrez_id",
    "gene_symbol",
    "id_repeats",
    "protein_id",
    "ptm_id",
    "ptm_peptide",
    "ratio_values",
    "ri_intensity",
    "tmt11_channel",
    "vial_label",
    "all_vial_labels",
    "tmt_plex"))

#' @title check proteomics ratio file
#'
#' @description check whether the proteomics ratio results files is following guidelines
#' @param df_ratio (data.frame) proteomics ratio results data frame (required)
#' @param isPTM (logical) `TRUE` if a PTM results file
#' @param f_proof (logical) `TRUE` (default) to print out distribution and NA plots
#' @param output_prefix (char) if `f_proof = TRUE`, provide a prefix for the output name
#' @param out_qc_folder (char) if `f_proof is TRUE`, a folder path can be provided
#' (otherwise print in current working directory)
#' @param printPDF (logical) if `TRUE` (default print plots to pdf)
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' test <- check_ratio_proteomics(df_ratio = metadata_metabolites_named,
#' isPTM =  TRUE, return_n_issues = TRUE, verbose = FALSE)
#' # "test" should be NULL
#' }
#' @export
check_ratio_proteomics <- function(df_ratio,
                                   isPTM,
                                   f_proof = TRUE,
                                   output_prefix = "ratio-file",
                                   out_qc_folder = NULL,
                                   printPDF = TRUE,
                                   return_n_issues = FALSE,
                                   verbose = TRUE){

  vial_label = all_vial_labels = tmt_plex = NULL

  ic_rr <- NULL

  if(isPTM){
    required_columns <- c("ptm_id", "protein_id", "gene_symbol", "entrez_id")
  }else{
    required_columns <- c("protein_id", "gene_symbol", "entrez_id")
  }

  if(!all( required_columns %in% colnames(df_ratio) )){
    if(verbose) message("      - (-) The following required columns are missed: ", appendLF = FALSE)
    if(verbose) message(paste(required_columns[!required_columns %in% colnames(df_ratio)], collapse = ", "))
    return(ic_rr)
  }

  ic_rr <- 0

  if(isPTM){
    # The ptm_id must be unique
    if( any(duplicated(df_ratio$ptm_id)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) NON UNIQUE ptm_id values ")
      ic <- ic + 1
      if(f_proof){
        # Print out redundancies
        # red_ids <- df_ratio$ptm_id[duplicated(df_ratio$ptm_id)]
        # rii_redundant_ids <- df_ratio[which(df_ratio$ptm_id %in% red_ids),]

        # Plot redundancies
        ya <- df_ratio %>% group_by(ptm_id) %>% summarise(id_repeats = n()) %>%
          group_by(id_repeats) %>% summarise(n = n()) %>%
          ggplot2::ggplot(aes(x = as.factor(id_repeats), y = n, fill = as.factor(id_repeats))) +
          geom_bar(stat = "identity") +
          theme_linedraw() +
          theme(legend.title = element_blank()) +
          labs(title = "RATIO RESULTS: Number of redundant ptm_ids",
               subtitle = paste(output_prefix))

        if(is.null(out_qc_folder)){
          out_plot_redundancies <- paste0(output_prefix,"-qc-ratior-redundant-ids.pdf")
        }else{
          out_plot_redundancies <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-ratior-redundant-ids.pdf"))
        }

        if(printPDF) pdf(out_plot_redundancies, width = 10, height = 8)
        print(ya)
        if(printPDF) garbage <- dev.off()
      }
    }else{
      if(verbose) message("   + (+) ptm_id unique: ok")
    }

    if( any(is.na(df_ratio$ptm_id)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) Some PTM_ID are missed")
    }else{
      if(verbose) message("   + (+) All PTM_ID available")
    }
  }else{
    if( any(duplicated(df_ratio$protein_id)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) Duplicated PROTEIN_ID found: FAIL")
      ic <- ic + 1
    }
  }

  if( any(is.na(df_ratio$protein_id)) ){
    ic_rr <- ic_rr + 1
    if(verbose) message("      - (-) Some PROTEIN_ID are missed")
  }else{
    if(verbose) message("   + (+) All PROTEIN_ID available")
  }

  if( any(is.na(df_ratio$gene_symbol)) ){
    ic_rr <- ic_rr + 1
    if(verbose) message("      - (-) Some GENE_SYMBOL are missed")
  }else{
    if(verbose) message("   + (+) All GENE_SYMBOL available")
  }

  if( any(is.na(df_ratio$entrez_id)) ){
    ic_rr <- ic_rr + 1
    if(verbose) message("      - (-) Some ENTREZ_ID are missed")
  }else{
    if(verbose) message("   + (+) All ENTREZ_ID available")
  }

  if(is.null(ic_rr) ){
    ic_rr <- 0
  }

  if(return_n_issues) return(ic_rr)
}

#' @title check proteomics reported ion intensity file
#'
#' @description check whether the proteomics rri results files is following guidelines
#' @param df_rri (data.frame) proteomics rri data frame (required)
#' @param isPTM (logical) `TRUE` if a PTM results file
#' @param f_proof (logical) `TRUE` (default) to print out distribution and NA plots
#' @param output_prefix (char) if `f_proof = TRUE`, provide a prefix for the output name
#' @param out_qc_folder (char) if `f_proof is TRUE`, a folder path can be provided
#' (otherwise print in current working directory)
#' @param printPDF (logical) if `TRUE` (default print plots to pdf)
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' test <- check_rii_proteomics(df_rri = metadata_metabolites_named,
#' isPTM =  TRUE, return_n_issues = TRUE, verbose = FALSE)
#' # "test" should be NULL
#' }
#' @export
check_rii_proteomics <- function(df_rri,
                                 isPTM,
                                 f_proof = TRUE,
                                 output_prefix = "rii-file",
                                 out_qc_folder = NULL,
                                 return_n_issues = FALSE,
                                 printPDF = TRUE,
                                 verbose = TRUE){


  ic_rii <- NULL

  if(isPTM){
    required_columns <- c("protein_id", "sequence", "ptm_id", "ptm_peptide", "gene_symbol", "entrez_id", "confident_score", "confident_site")
  }else{
    required_columns <- c("protein_id", "sequence", "gene_symbol", "entrez_id")
  }

  if(!all( required_columns %in% colnames(df_rri) )){
    if(verbose) message("      - (-) The following required columns are missed: ", appendLF = FALSE)
    if(verbose) message(paste(required_columns[!required_columns %in% colnames(df_rri)], collapse = ", "))
    return(ic_rii)
  }

  ic_rii <- 0

  if(isPTM){
    # The ptm_id must be unique
    if( any(duplicated(df_rri$ptm_peptide)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) NON UNIQUE PTM_PEPTIDE values")
      ic <- ic + 1
      if(f_proof){
        # Print out redundancies
        # red_ids <- df_rri$ptm_peptide[duplicated(df_rri$ptm_peptide)]
        # rii_redundant_ids <- df_rri[which(df_rri$ptm_peptide %in% red_ids),]

        # Plot redundancies
        df <- df_rri %>%
          dplyr::group_by(ptm_peptide) %>%
          dplyr::summarise(id_repeats = n()) %>%
          dplyr::group_by(id_repeats) %>%
          dplyr::summarise(n = n())

        xa <- ggplot2::ggplot(df, aes(x = as.factor(id_repeats), y = n, fill = as.factor(id_repeats))) +
          geom_bar(stat = "identity") +
          theme_linedraw() +
          theme(legend.title = element_blank()) +
          labs(title = "RII: Number of redundant ptm_peptides",
               subtitle = paste(output_prefix))
        if(verbose) message("      - (p) Plot redundant values")

        if(is.null(out_qc_folder)){
          out_plot_redundancies <- paste0(output_prefix,"-qc-redundant-ids.pdf")
        }else{
          out_plot_redundancies <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-redundant-ids.pdf"))
        }

        if(printPDF) pdf(out_plot_redundancies, width = 10, height = 8)
        print(xa)
        if(printPDF) garbage <- dev.off()
      }
    }

    if( any(is.na(df_rri$sequence)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) Some SEQUENCE are missed")
    }else{
      if(verbose) message("   + (+) All SEQUENCE available")
    }

    if( any(is.na(df_rri$ptm_id)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) Some PTM_ID are missed")
    }else{
      if(verbose) message("   + (+) All PTM_ID available")
    }

    if( any(is.na(df_rri$confident_score)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) Some CONFIDENT_SCORE are missed")
    }else{
      if(verbose) message("   + (+) All CONFIDENT_SCORE available")
    }

    if( any(is.na(df_rri$confident_site)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) Some CONFIDENT_SITE are missed")
    }else{
      if(verbose) message("   + (+) All CONFIDENT_SITE available")
    }
  }else{
    df_rri$protein_sequence <- paste0(df_rri$protein_id,"-", df_rri$sequence)
    if(any(duplicated(df_rri$protein_sequence))){
      ic_rii <- ic_rii + 1
      ic <- ic + 1
      if(verbose) message("      - (-) Duplicated Protein + Sequence identified")
    }else{
      if(verbose) message("   + (+) All Protein + Sequence are unique")
    }
    df_rri$protein_sequence <- NULL
  }

  if( any(is.na(df_rri$protein_id)) ){
    ic_rii <- ic_rii + 1
    if(verbose) message("      - (-) Some PROTEIN_ID are missed")
  }else{
    if(verbose) message("   + (+) All PROTEIN_ID available")
  }

  if( any(is.na(df_rri$gene_symbol)) ){
    ic_rii <- ic_rii + 1
    if(verbose) message("      - (-) Some GENE_SYMBOL are missed")
  }else{
    if(verbose) message("   + (+) All GENE_SYMBOL available")
  }

  if( any(is.na(df_rri$entrez_id)) ){
    ic_rii <- ic_rii + 1
    if(verbose) message("      - (-) Some ENTREZ_ID are missed")
  }else{
    if(verbose) message("   + (+) All ENTREZ_ID available")
  }

  if( is.null(ic_rii) ){
    ic_rii <- 0
  }

  if(return_n_issues) return(ic_rii)
}


#' @title check proteomics vial metadata file
#'
#' @description check whether the proteomics rri results files is following guidelines
#' @param df_vm (data.frame) proteomics vial_label data frame (required)
#' @param return_n_issues (logical) if `TRUE` returns the number of issues.
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' test <- check_vial_metadata_proteomics(df_vm = metadata_metabolites_named,
#' return_n_issues = TRUE, verbose = FALSE)
#' # "test" should be NULL
#' }
#' @export
check_vial_metadata_proteomics <- function(df_vm,
                                           return_n_issues = FALSE,
                                           verbose = TRUE){

  ic_vm <- NULL

  required_columns <- c("vial_label", "tmt_plex", "tmt11_channel")

  if(!all( required_columns %in% colnames(df_vm) )){
    if(verbose) message("      - (-) The following required columns are missed: ", appendLF = FALSE)
    if(verbose) message(paste(required_columns[!required_columns %in% colnames(df_vm)], collapse = ", "))
    return(ic_vm)
  }

  ic_vm <- 0

  df_vm$vial_label <- gsub(" ", "", df_vm$vial_label)
  valid_channels <- c("126C", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C")
  plexes <- unique(df_vm$tmt_plex)

  for(p in plexes){
    temp_plex <- df_vm[which(df_vm$tmt_plex == p),]
    if(all(valid_channels %in% temp_plex$tmt11_channel)){
      if(verbose) message("   + (+) All tmt channels are valid in plex ", paste(p))
    }else{
      if(verbose) message("      - (-) Invalid tmt channels in the file")
      ic_vm <- ic_vm + 1
    }
  }

  all_samples <- df_vm$vial_label
  all_vial_labels <- NA

  if( any( grepl("Ref", df_vm$vial_label) ) ){
    all_vial_labels <- all_samples[!grepl('^Ref', all_samples)]
  }else{
    if(verbose) message("      - (-) Ref channels not found in vial_metadata")
    ic_vm <- ic_vm + 1
  }

  if(is.null(ic_vm)){
    ic_vm <- 0
  }

  if(return_n_issues) return(ic_vm)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Load Proteomics batch
#'
#' @description Load a proteomics batch
#' @param input_results_folder (char) path to the PROCESSED folder to check
#' @param isPTM (logical) `TRUE` if it is Post-Translational-Modification proteomics assay
#' @param verbose (logical) `TRUE` (default) prints QC details.
#' @return (List) of data frames: rii,
#' @export
load_proteomics <- function(input_results_folder,
                            isPTM,
                            verbose = TRUE){


  if(any(missing(input_results_folder) |
         missing(isPTM)))
    stop("One (or many) of the required arguments missed.
        Please, check the help for this function to find out more")

  # Validate folder structure-----
  processfolder <- validate_processFolder(input_results_folder)
  phase <- validate_phase(input_results_folder)

  if(verbose) message("# LOAD PROTEOMICS BATCH\n\n")
  if(verbose) message(paste(processfolder))

  # COUNTS----
  ic <- 0 # Critical issues
  ic_rii <- NA # RII issues
  ic_rr <- NA # ratio results
  ic_vm <- NA # vial metadata
  ic_man <- NA # namifiest


  # VIAL METADATA-----
  if(verbose) message("\n## VIAL METADATA\n")

  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "vial_metadata.txt",
                     verbose = verbose)
  f_vm <- lista$flag

  all_vial_labels <- NA

  if(f_vm){
    v_m <- lista$df
    v_m <- filter_required_columns(df = v_m,
                                   type = "v_m",
                                   verbose = verbose)
    # Remove space between "Ref_ A"
    v_m$vial_label <- gsub(" ", "", v_m$vial_label)
    ic_vm <- check_vial_metadata_proteomics(df_vm = v_m,
                                            return_n_issues = TRUE,
                                            verbose = verbose)

    if(is.null(ic_vm)){
      f_vm <- FALSE
      ic <- ic + 1
      ic_vm <- NA
    }else{
      if( any( grepl("Ref", v_m$vial_label) ) ){
        all_samples <- v_m$vial_label
        all_vial_labels <- all_samples[!grepl('^Ref', all_samples)]
      }else{
        if(verbose) message("      - (-) Ref channels not found in vial_metadata")
        ic_vm <- ic_vm + 1
      }
    }
  }else{
    if(verbose) message("      - (-) {vial_metadata} not available")
    ic <- ic + 1
  }


  # REPORTED ION INTENSITY -----
  if(verbose) message("\n## REPORTED ION INTENSITY\n")

  if(verbose) message("   + Loading the file (might take some time)")
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results_RII-peptide.txt",
                     verbose = verbose)
  f_rii <- lista$flag
  if(f_rii){
    peprii <- lista$df
    ic_rii <- check_rii_proteomics(df_rri = peprii,
                                   isPTM = isPTM,
                                   verbose = verbose)

  }else{
    if(verbose) message("      - (-) {results_RII-peptide} file not available")
    ic <- ic + 1
  }

  # RATIO----

  if(verbose) message("\n\n## RATIO results\n")

  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results_ratio.txt",
                     verbose = verbose)

  f_rr <- lista$flag
  if(f_rr){
    ratior <- lista$df
    ic_rr <- check_ratio_proteomics(df_ratio = ratior,
                                    isPTM = isPTM,
                                    verbose = verbose)
  }else{
    if(verbose) message("      - (-) {results_ratio.txt} file not available")
    ic <- ic + 1
  }

  if(f_rr & f_rii & f_vm){
    list_df <- list ("peprii" = peprii,
                     "v_m" = v_m,
                     "ratior" = ratior,
                     "phase" = phase)
    return(list_df)
  }else{
    stop("PROBLEMS LOADING PROTEOMICS BATCH")
  }

}







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Validate a Proteomics submissions
#'
#' @description Validate a Proteomics submission
#' @param input_results_folder (char) path to the PROCESSED folder to check
#' @param isPTM (logical) `TRUE` if it is Post-Translational-Modification proteomics assay
#' @param cas (char) CAS code
#' @param dmaqc_shipping_info (char) phase code
#' @param f_proof (char) print out pdf with charts including:
#' - Reported Ion Intensity boxplot distribution and percentage of NA values per sample
#' - Ratio: ratio boxplot distribution and percentage of NA values per samples
#' @param out_qc_folder (char) if f_proof is TRUE, then a folder must be provided
#' @param return_n_issues (logical) if `TRUE` (default) returns the number of issues
#' @param full_report (logical) if `FALSE` (default) it returns only the
#' total number of issues. If `TRUE` returns the details of the number of issues (by
#' group of files, e.g., results, metadata_metabolites, etc)
#' @param printPDF (logical) if `TRUE` (default print plots to pdf)
#' @param verbose (logical) `TRUE` (default) prints QC details.
#' @return (data.frame) Summary of issues
#' @export
validate_proteomics <- function(input_results_folder,
                                isPTM,
                                cas,
                                dmaqc_shipping_info = NULL,
                                f_proof = FALSE,
                                out_qc_folder = NULL,
                                return_n_issues = TRUE,
                                full_report = FALSE,
                                printPDF = TRUE,
                                verbose = TRUE){

  if(any(missing(input_results_folder) |
         missing(isPTM) |
         missing(cas)))
    stop("One (or many) of the required arguments missed.
        Please, check the help for this function to find out more")

  # Validate folder structure-----
  processfolder <- validate_processFolder(input_results_folder)
  assay <- validate_assay(input_results_folder)
  phase <- validate_phase(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)

  # Print out proofs----
  if(f_proof){
    if(!is.null(out_qc_folder)){
      if(!dir.exists(file.path(out_qc_folder))){
        dir.create(file.path(out_qc_folder), recursive = TRUE)
      }
    }

    output_prefix <- paste0(cas, ".", tolower(phase), ".", tissue_code, ".",tolower(assay), ".", tolower(processfolder))
  }

  input_folder_short <- regmatches(input_results_folder, regexpr("PASS.*RESULTS_[0-9]{8}", input_results_folder))
  if(purrr::is_empty(input_folder_short)){
    if(verbose) message("\nThe RESULTS_YYYYMMDD folder full path is not correct. Example:")
    if(verbose) message("/full/path/to/folder/PASS1A-06/T66/RPNEG/BATCH1_20190822/RESULTS_202003")
    stop("Input folder not according to guidelines")
  }

  if(verbose) message("# PROTEOMICS QC report\n\n")
  if(verbose) message("+ Site: ", toupper(cas))
  if(verbose) message("+ Folder: `",paste0(input_folder_short),"`")

  # COUNTS----
  ic <- 0 # Critical issues
  ic_rii <- NA # RII issues
  ic_rr <- NA # ratio results
  ic_vm <- NA # vial metadata
  ic_man <- NA # manifest

  if(missing(dmaqc_shipping_info)){
    dmaqc_shipping_info = NULL
  }

  if(is.null(dmaqc_shipping_info)){
    ic_vl <- "missed"
  } else{
    ic_vl <- NA
  }

  if(missing(out_qc_folder)){
    out_qc_folder = NULL
  }

  # if(is.null(dmaqc_shipping_info)){
  #   message("dmaqc_shipping_info is ", paste(dmaqc_shipping_info), " imagino que NULL\n\n")
  # }
  #
  # if(is.null(out_qc_folder)){
  #   message("out_qc_folder is ", paste(out_qc_folder), " imagino que NULL\n\n")
  # }
  #

  #
  # message("ic_vl is now ", paste(ic_vl), " or ", ic_vl)

  # VIAL METADATA-----

  if(verbose) message("\n## VIAL METADATA\n")

  lista <- open_file(input_results_folder = input_results_folder,
                                   filepattern = "vial_metadata.txt",
                                   verbose = verbose)
  f_vm <- lista$flag

  all_vial_labels <- NA

  if(f_vm){
    v_m <- lista$df
    v_m <- filter_required_columns(df = v_m,
                                   type = "v_m",
                                   verbose = verbose)
    ic_vm <- check_vial_metadata_proteomics(df_vm = v_m,
                                            return_n_issues = TRUE,
                                            verbose = verbose)
    if(is.null(ic_vm)){
      f_vm <- FALSE
      ic <- ic + 1
      ic_vm <- NA
    }else{
      if( any( grepl("Ref", v_m$vial_label) ) ){
        all_samples <- v_m$vial_label
        all_vial_labels <- all_samples[!grepl('^Ref', all_samples)]
      }else{
        if(verbose) message("      - (-) Ref channels not found in vial_metadata")
        ic_vm <- ic_vm + 1
      }
    }
  }else{
    if(verbose) message("      - (-) {vial_metadata} not available")
    ic <- ic + 1
  }


  # REPORTED ION INTENSITY -----
  if(verbose) message("\n## REPORTED ION INTENSITY\n")

  if(verbose) message("   + Loading the file (might take some time)")
  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results_RII-peptide.txt",
                     verbose = verbose)

  f_rii <- lista$flag
  if(f_rii){
    peprii <- lista$df
    ic_rii <- check_rii_proteomics(df_rri = peprii,
                                  isPTM = isPTM,
                                  f_proof = f_proof,
                                  output_prefix = output_prefix,
                                  return_n_issues = TRUE,
                                  out_qc_folder = out_qc_folder,
                                  printPDF = printPDF,
                                  verbose = verbose)

    if(is.null(ic_rii)){
      f_rii <- FALSE
      ic_rii <- NA
      ic <- ic + 1
    }else{
      if(f_proof){

        if(verbose) message("   + (+) PLOT: RII distribution and NA values")

        # Check distributions
        if(isPTM){
          peptides_long <- peprii %>% tidyr::pivot_longer(cols = -c(ptm_peptide, ptm_id, protein_id, gene_symbol, sequence, entrez_id),
                                                          names_to = "vial_label",
                                                          values_to = "ri_intensity")
        }else{
          peptides_long <- peprii %>% tidyr::pivot_longer(cols = -c(protein_id, gene_symbol, sequence, entrez_id),
                                                          names_to = "vial_label",
                                                          values_to = "ri_intensity")
        }

        # Only if vial_label metadata is available
        if(f_vm){

          if(verbose) message("   + (p) Plot intensity distributions")

          peptides_long <- merge(v_m, peptides_long, by = c("vial_label"))

          peptides_long$vial_label <- as.character(peptides_long$vial_label)
          peptides_long$vial_label <- as.factor(peptides_long$vial_label)
          peptides_long$tmt_plex <- as.factor(peptides_long$tmt_plex)

          pise <- ggplot2::ggplot(peptides_long,
                                  aes(x = reorder(vial_label, log2(ri_intensity), FUN = median, na.rm = TRUE),
                                      y = log2(ri_intensity),
                                      fill = tmt_plex)) +
            geom_boxplot(na.rm = TRUE) +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90,
                                             hjust = 1,
                                             vjust = 0.5,
                                             size = 8)) + #legend.position = "none"
            labs(x = "vial_label", y = "log2(ri_intensity)") +
            labs(title = "Reporter Ion Intensity",
                 subtitle = output_prefix)


          # Plottingh NA values
          if(verbose) message("   + (p) Plot NA values")

          p_na_peprii <- peprii %>%
            inspectdf::inspect_na() %>%
            dplyr::arrange(match(col_name, colnames(peprii))) %>%
            inspectdf::show_plot() +
            ylim(0, 100) + theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90,
                                             hjust = 1,
                                             vjust = 0.5,
                                             size = 8))

          if(is.null(out_qc_folder)){
            out_plot_dist <- paste0(output_prefix,"-qc-rii-distribution.pdf")
          }else{
            out_plot_dist <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-rii-distribution.pdf"))
          }

          if(printPDF) pdf(out_plot_dist, width = 12, height = 8)
          print(pise)
          print(p_na_peprii)
          if(printPDF) garbage <- dev.off()
        }
      } # print plots
    }
  }else{
    if(verbose) message("      - (-) {results_RII-peptide} file not available")
    ic <- ic + 1
  }

  # RATIO----

  if(verbose) message("\n\n## RATIO results\n")

  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "results_ratio.txt",
                     verbose = verbose)

  f_rr <- lista$flag
  if(f_rr){
    ratior <- lista$df
    ic_rr <- check_ratio_proteomics(df_ratio = ratior,
                                    isPTM = isPTM,
                                    f_proof = f_proof,
                                    output_prefix = output_prefix,
                                    out_qc_folder = out_qc_folder,
                                    printPDF = printPDF,
                                    return_n_issues = TRUE,
                                    verbose = verbose)

    if(is.null(ic_rr)){
      f_rr <- FALSE
      ic <- ic + 1
      ic_rr <- NA
    }else{
      # Plot distributions

      if(f_proof){

        if(verbose) message("   + (+) PLOT: RATIO distribution and NA values")

        # Check distributions
        if(isPTM){
          ratior_long <- ratior %>% tidyr::pivot_longer(cols = -c(ptm_id, protein_id, gene_symbol, entrez_id, confident_score, confident_site),
                                                            names_to = "vial_label",
                                                            values_to = "ratio_values")
        }else{
          ratior_long <- ratior %>% tidyr::pivot_longer(cols = -c(protein_id, gene_symbol, entrez_id),
                                                            names_to = "vial_label",
                                                            values_to = "ratio_values")
        }

        if(f_vm){
          if(verbose) message("   + (p) Plotting ratio distributions")

          ratior_long <- merge(v_m, ratior_long, by = c("vial_label"))

          pisr <- ggplot2::ggplot(ratior_long,
                                  aes(x = reorder(vial_label, ratio_values, FUN = median, na.rm = TRUE),
                                      y = ratio_values,
                                      fill = tmt_plex)) +
            geom_boxplot(na.rm = TRUE) +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90,
                                             hjust = 1,
                                             vjust = 0.5,
                                             size = 8)) + #legend.position = "none"
            labs(x = "vial_label", y = "ratio_values") +
            labs(title = "Ratio",
                 subtitle = output_prefix)

          # Plotting NA values
          if(verbose) message("   + (p) Plotting NA percentage in ratio results")
          p_na_ratior <- ratior %>%
            inspectdf::inspect_na() %>%
            dplyr::arrange(match(col_name, colnames(ratior))) %>%
            inspectdf::show_plot() +
            ylim(0, 100) + theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90,
                                             hjust = 1,
                                             vjust = 0.5,
                                             size = 8))

          if(is.null(out_qc_folder)){
            out_plot_ratdist <- paste0(output_prefix,"-qc-ratio-distribution.pdf")
          }else{
            out_plot_ratdist <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-ratio-distribution.pdf"))
          }

          if(printPDF) pdf(out_plot_ratdist, width = 12, height = 8)
          print(pisr)
          print(p_na_ratior)
          if(printPDF) garbage <- dev.off()
        }
      }
    } #print plots
  }else{
    if(verbose) message("      - (-) {results_ratio.txt} file not available")
    ic <- ic + 1
  }

  # MANIFEST----

  if(verbose) message("\n## MANIFEST\n")

  batch <- gsub("(.*)(RESULTS.*)", "\\1", input_results_folder)

  file_manifest <- list.files(normalizePath(batch),
                                     pattern="file_manifest",
                                     ignore.case = TRUE,
                                     full.names=TRUE,
                                     recursive = TRUE)

  if(length(file_manifest) == 0){
    f_man <- FALSE
    ic_man <- 1
  }else if(length(file_manifest) > 1){
    file_manifest <- file_manifest[length(file_manifest)]
    f_man <- TRUE
  }else if( length(file_manifest) == 1 ){
    f_man <- TRUE
  }

  if(f_man){
    manifest <- read.csv(file_manifest, stringsAsFactors = FALSE)
    mani_columns <- c("file_name", "md5")
    if( all(colnames(manifest) %in% mani_columns ) ){
      if(verbose) message("   + (+)  (file_name, md5) columns available in manifest file")
      if(f_rr){
        ratio_file <- manifest$file_name[grepl("ratio", manifest$file_name)]
        if(file.exists(file.path(batch, ratio_file))){
          if(verbose) message("   + (+) ratio file included")
        }else{
          if(verbose) message("      - (-) ratio file is not included in manifest file")
          ic_man <- 1
        }
      }

      if(f_rii){
        rii_file <- manifest$file_name[grepl("RII", manifest$file_name)]
        if(file.exists(file.path(batch, rii_file))){
          if(verbose) message("   + (+) RII file included")
        }else{
          if(verbose) message("      - (-) RII file is not included in manifest file")
          ic_man <- 1
        }
      }

      if(f_vm){
        vm_file <- manifest$file_name[grepl("vial_metadata", manifest$file_name)]
        if(file.exists(file.path(batch, vm_file))){
          if(verbose) message("   + (+) VIAL_METADATA file included")
        }else{
          if(verbose) message("      - (-) VIAL_METADATA file is not included in manifest file")
          ic_man <- 1
        }
      }

      if( any(is.na(manifest$md5)) ){
        if(verbose) message("      - (-) MD5 column contains NA values")
      }

    }else{
      if(verbose) message("      - (-) Not all the columns are available")
      ic_man <- ic_man + 1
      ic <- ic + 1
    }
  }

  # CHECK DMAQC----

  # Validate vial labels from DMAQC

  if(verbose) message("\n\n## DMAQC validation\n")
  failed_samples <- check_failedsamples(input_results_folder = input_results_folder, verbose = verbose)

  if( is.na(ic_vl) ){
    if(f_vm){
      message("before check viallabel, tell me the valud of ic_vl: ", paste(ic_vl),"\n\n")
      ic_vl <- check_viallabel_dmaqc(vl_submitted = all_vial_labels,
                                     tissue_code = tissue_code,
                                     cas = cas,
                                     phase = phase,
                                     failed_samples = failed_samples,
                                     dmaqc_shipping_info = dmaqc_shipping_info,
                                     return_n_issues = TRUE,
                                     verbose = verbose)
      message("After")
    }else{
      ic <- ic + 1
      if(verbose) message("      - (-) DMAQC validation cannot be performed (no vial label data available)")
    }
  }

  # INTER-FILE VALIDATION----

  if(verbose) message("\n## INTER-file validations\n")

  if(f_vm & f_rii & f_rr){

    if( all(all_vial_labels %in% colnames(peprii)) ){
      if(verbose) message("   + (+) All all_vial_labels in RII file")
    }else{
      ic <- ic + 1
      if(verbose) message("      - (-) Some all_vial_labels not available in RII file")
    }

    if( all(all_vial_labels %in% colnames(ratior)) ){
      if(verbose) message("   + (+) All all_vial_labels in RATIO results file")
    }else{
      ic <- ic + 1
      if(verbose) message("      - (-) Some all_vial_labels not available in RATIO results file")
    }

  }else{
    if(verbose) message("      - (-) One of the files is not available")
  }

  # PRINT OUT RESULTS-----
  batchversion <- stringr::str_extract(string = input_results_folder, pattern = "BATCH.*_[0-9]+/RESULTS_[0-9]+")

  qc_date <- Sys.time()
  qc_date <- gsub("-", "", qc_date)
  qc_date <- gsub(" ", "_", qc_date)
  qc_date <- gsub(":", "", qc_date)
  t_name <- bic_animal_tissue_code$bic_tissue_name[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)]

  if(return_n_issues){
    total_issues <- sum(ic, ic_vm, ic_rr, ic_rii, na.rm = TRUE)
    if(verbose) message("\nTOTAL NUMBER OF ISSUES: ", total_issues,"\n")
    if(full_report){
      reports <- data.frame(cas = cas,
                            phase= phase,
                            tissue = tissue_code,
                            t_name = t_name,
                            assay = assay,
                            version = batchversion,
                            vial_label = length(all_vial_labels),
                            qc_samples = NA,
                            dmaqc_valid = ic_vl,
                            critical_issues = ic,
                            pep_rri = ic_rii,
                            ratio = ic_rr,
                            vial_meta = ic_vm,
                            qc_date = qc_date)
      return(reports)
    }else{
      return(total_issues)
    }
  }
}



#' @title Write proteomics data release
#'
#' @description Write out proteomics data releases. Doesn't check whether
#' data has been submited according to guidelines
#' @param input_results_folder (char) Path to the PROCESSED_YYYYMMDD folder
#' @param isPTM (logical) is a proteomics PTM experiment?
#' @param folder_name (char) output folder name.
#' @param folder_root (char) absolute path to write the output folder. Default: current directory
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return bic release folder/file structure `PHASE/OMICS/TCODE_NAME/ASSAY/` and file names, including:
#' `motrpac_YYYYMMDD_phasecode_tissuecode_omics_assay_file-details.txt` where files-details can be:
#' `rii-results.txt`, `ratio-results.txt`, `vial-metadata.txt`
#' @examples
#' \dontrun{
#' write_proteomics_releases(
#'    input_results_folder = "/full/path/to/RESULTS_YYYYMMDD/")
#' }
#' @export
write_proteomics_releases <- function(input_results_folder,
                                      isPTM,
                                      folder_name = "motrpac_release",
                                      folder_root = NULL,
                                      verbose = TRUE){


  # Get names from input_results_folder------
  assay <- validate_assay(input_results_folder)
  phase <- validate_phase(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)
  folder_phase <- tolower(phase)
  folder_tissue <- bic_animal_tissue_code$tissue_name_release[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)]

  assay_codes$cas_code <- NULL
  assay_codes <- unique(assay_codes)
  if(length(assay_codes$assay_code[which(assay_codes$submission_code == assay)]) == 1){
    folder_assay <- assay_codes$assay_code[which(assay_codes$submission_code == assay)]
  }else{
    stop("ASSAY code ", assay, " not available in < assay_codes >")
  }

  if(verbose) message("+ Writing out ", phase, " ", tissue_code, " ", assay, " files", appendLF = FALSE)

  # Load all datasets----
  prot_dfs <- load_proteomics(input_results_folder = input_results_folder,
                              isPTM = isPTM,
                              verbose = FALSE)

  # Create output folder-------
  if (is.null(folder_root)){
    folder_root <- getwd()
  }else{
    folder_root <- normalizePath(folder_root)
  }

  output_folder <- file.path(folder_root, folder_name, folder_phase, "proteomics", folder_tissue, folder_assay)

  if(!dir.exists(file.path(output_folder))){
    dir.create(file.path(output_folder), recursive = TRUE)
  }

  file_name_shared <- paste0("motrpac_",
                             folder_phase, "_",
                             folder_tissue, "_",
                             folder_assay)


  # Create and write FILES-----
  prot_rii <- file.path(output_folder, paste0(file_name_shared,"_rii-results.txt"))
  prot_ratio <- file.path(output_folder, paste0(file_name_shared,"_ratio-results.txt"))
  vial_metadata <- file.path(output_folder, paste0(file_name_shared,"_vial-metadata.txt"))

  write.table(prot_dfs$peprii, prot_rii, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(prot_dfs$ratior, prot_ratio, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(prot_dfs$v_m, vial_metadata, row.names = FALSE, sep = "\t", quote = FALSE)

  if(verbose) message("...done!")
}

