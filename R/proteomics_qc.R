
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

  ic_rr <- 0

  required_columns <- get_required_columns(isPTM = isPTM,
                                           prot_file = "ratio")

  if(!all( required_columns %in% colnames(df_ratio) )){
    if(verbose) message("      - (-) CRITICAL REQUIRED column(s) missed: ", appendLF = FALSE)
    if(verbose) message(paste(required_columns[!required_columns %in% colnames(df_ratio)], collapse = ", "))
    ic_rr <- 10
    return(ic_rr)
  }

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

    if( any(is.na(df_ratio$confident_score)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) Some CONFIDENT_SCORE are missed")
    }else{
      if(verbose) message("   + (+) All CONFIDENT_SCORE available")
    }

    if( any(is.na(df_ratio$confident_site)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) Some CONFIDENT_SITE are missed")
    }else{
      if(verbose) message("   + (+) All CONFIDENT_SITE available")
    }
    
    if( any(is.na(df_ratio$flanking_sequence)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) Some FLANKING_SEQUENCEs are missed")
    }else{
      if(verbose) message("   + (+) All FLANKING_SEQUENCEs available")
    }
    
    if("ptm_score" %in% colnames(df_ratio)){
      if( any(is.na(df_ratio$ptm_score)) ){
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) Some PTM_SCOREs are missed")
      }else{
        if(verbose) message("   + (+) All PTM_SCOREs available")
      }
      if(is.numeric(df_ratio$ptm_score)){
        if(verbose) message("   + (+) PTM_SCORE is numeric")
      }else{
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) PTM_SCORE is  NOT NUMERIC")
      }
    }else{
      ic_rr <- ic_rr + 2
      if(verbose) message("      - (-) ptm_score column missed")
    }
  }else{
    if( any(duplicated(df_ratio$protein_id)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) Duplicated PROTEIN_ID found: FAIL")
      ic <- ic + 1
    }
    
    if("num_peptides" %in% colnames(df_ratio)){
      if( any(is.na(df_ratio$num_peptides)) ){
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) Some num_peptides are missed")
      }else{
        if(verbose) message("   + (+) All num_peptides available")
      }
      if(is.numeric(df_ratio$num_peptides)){
        if(verbose) message("   + (+) num_peptides is numeric")
      }else{
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) num_peptides is NOT NUMERIC")
      }
    }else{
      ic_rr <- ic_rr + 2
      if(verbose) message("      - (-) num_peptides column missed")
    }
    
    if("percent_coverage" %in% colnames(df_ratio)){
      if( any(is.na(df_ratio$percent_coverage)) ){
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) Some percent_coverage are missed")
      }else{
        if(verbose) message("   + (+) All percent_coverage available")
      }
      if(is.numeric(df_ratio$percent_coverage)){
        if(verbose) message("   + (+) percent_coverage is numeric")
      }else{
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) percent_coverage is NOT NUMERIC")
      }
    }else{
      ic_rr <- ic_rr + 2
      if(verbose) message("      - (-) percent_coverage column missed")
    }
    
    if("protein_score" %in% colnames(df_ratio)){
      if( any(is.na(df_ratio$protein_score)) ){
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) Some protein_score are missed")
      }else{
        if(verbose) message("   + (+) All protein_score available")
      }
      if(is.numeric(df_ratio$protein_score)){
        if(verbose) message("   + (+) protein_score is numeric")
      }else{
        ic_rr <- ic_rr + 1
        if(verbose) message("      - (-) protein_score is NOT NUMERIC")
      }
    }else{
      ic_rr <- ic_rr + 2
      if(verbose) message("      - (-) protein_score column missed")
    }

  } # PROTEIN RATIO

  if( any(is.na(df_ratio$protein_id)) ){
    ic_rr <- ic_rr + 1
    if(verbose) message("      - (-) Some PROTEIN_ID are missed")
  }else{
    if(verbose) message("   + (+) All PROTEIN_ID available")
  }

  if( any( is.na( df_ratio$gene_symbol[which(df_ratio$is_contaminant == FALSE)] ) ) ){
    # ic_rii <- ic_rii + 1
    
    if(verbose) message(paste("      - (!) GENE_SYMBOL values are missed for", 
                              length( df_ratio$protein_id[which(is.na( df_ratio$gene_symbol[which(df_ratio$is_contaminant == FALSE)] ))] ), 
                              "protein ids"))
    
  }else{
    if(verbose) message("   + (+) All GENE_SYMBOL available")
  }

  if( any( is.na( df_ratio$entrez_id[which(df_ratio$is_contaminant == FALSE)] ) ) ){
    # ic_rii <- ic_rii + 1
    if(verbose) message(paste("      - (!) ENTREZ_ID values are missed for", 
                              length(df_ratio$protein_id[which(is.na( df_ratio$entrez_id[which(df_ratio$is_contaminant == FALSE)] ))]), 
                              "protein ids"))
  }else{
    if(verbose) message("   + (+) All ENTREZ_ID available")
  }
  
  if( any(is.na(df_ratio$redundant_ids)) ){
    ic_rr <- ic_rr + 1
    if(verbose) message("      - (-) Some redundant_ids are missed")
  }else{
    if(verbose) message("   + (+) All redundant_ids available")
  }
  
  if( any(is.na(df_ratio$organism_name)) ){
    ic_rr <- ic_rr + 1
    if(verbose) message("      - (-) Some organism_name are missed")
  }else{
    if(verbose) message("   + (+) All organism_name available")
  }
  
  if("is_contaminant" %in% colnames(df_ratio)){
    if(is.logical(df_ratio$is_contaminant)){
      if(verbose) message("   + (+) All is_contaminant available")
    }else{
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) is_contaminant is not logical (TRUE/FALSE)")
    }
    if( any(is.na(df_ratio$is_contaminant)) ){
      ic_rr <- ic_rr + 1
      if(verbose) message("      - (-) Some is_contaminant are missed")
    }else{
      if(verbose) message("   + (+) All is_contaminant available")
    }
  }
  
  if( any(is.na(df_ratio$is_contaminant)) ){
    ic_rr <- ic_rr + 1
    if(verbose) message("      - (-) Some organism_name are missed")
  }else{
    if(verbose) message("   + (+) All organism_name available")
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


  ic_rii <- 0

  required_columns <- get_required_columns(isPTM = isPTM,
                                           prot_file = "rii")

  if(!all( required_columns %in% colnames(df_rri) )){
    ic_rii <- 10
    if(verbose) message("      - (-) Error! Critical Required column(s) not found: ", appendLF = FALSE)
    if(verbose) message(paste(required_columns[!required_columns %in% colnames(df_rri)], collapse = ", "))
    return(ic_rii)
  }

  if(isPTM){
    # The ptm_id must be unique
    if( any(duplicated(df_rri$ptm_peptide)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) ERROR: NON UNIQUE PTM_PEPTIDE values")
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

    if( any(is.na(df_rri$ptm_id)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) ERROR: Some PTM_ID values are missed")
    }else{
      if(verbose) message("   + (+) All PTM_ID available")
    }
    
    if("confident_score" %in% colnames(df_rri)){
      if(is.numeric(df_rri$confident_score)){
        if(verbose) message("   + (+) All CONFIDENT_SCORE are numeric")
      }else{
        ic_rii <- ic_rii + 1
        if(verbose) message("      - (-) ERROR: CONFIDENT_SCORE is not numeric!")
      }
      
      if( any(is.na(df_rri$confident_score)) ){
        ic_rii <- ic_rii + 1
        if(verbose) message("      - (-) ERROR: Some CONFIDENT_SCORE values are missed")
      }else{
        if(verbose) message("   + (+) All CONFIDENT_SCORE available")
      }
    }

    if("confident_site" %in% colnames(df_rri)){
      if( any(is.na(df_rri$confident_site)) ){
        ic_rii <- ic_rii + 1
        if(verbose) message("      - (-) ERROR: Some CONFIDENT_SITE are missed")
      }else{
        if(verbose) message("   + (+) All CONFIDENT_SITE available")
      }
      
      if(is.logical(df_rri$confident_site)){
        if(verbose) message("   + (+) confident_site is a logical variable")
      }else{
        ic_rii <- ic_rii + 1
        if(verbose) message("   + (+) confident_site is NOT a logical variable (all values should be TRUE/FALSE)")
      }
    }else{
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) ERROR: confident_site column is not available")
    }
    # isPTM == TRUE
  }else{
    df_rri$protein_sequence <- paste0(df_rri$protein_id,"-", df_rri$sequence)
    if(any(duplicated(df_rri$protein_sequence))){
      ic_rii <- ic_rii + 1
      ic <- ic + 1
      if(verbose) message("      - (-) ERROR: Duplicated Protein + Sequence identified")
    }else{
      if(verbose) message("   + (+) All Protein + Sequence are unique")
    }
    df_rri$protein_sequence <- NULL
  }

  if( any(is.na(df_rri$protein_id)) ){
    ic_rii <- ic_rii + 1
    if(verbose) message("      - (-) ERROR: (some) PROTEIN_IDs missed")
  }else{
    if(verbose) message("   + (+) All PROTEIN_ID values available")
  }

  if( any( is.na( df_rri$gene_symbol[which(df_rri$is_contaminant == FALSE)] ) ) ){
    # ic_rii <- ic_rii + 1
    if(verbose) message(paste("      - (!) GENE_SYMBOL values are missed for", 
                              length(df_rri$protein_id[which(is.na( df_rri$gene_symbol[which(df_rri$is_contaminant == FALSE)] ))]), 
                              "protein ids"))
  }else{
    if(verbose) message("   + (+) All GENE_SYMBOL ids available")
  }

  if( any( is.na( df_rri$entrez_id[which(df_rri$is_contaminant == FALSE)] ) ) ){
    # ic_rii <- ic_rii + 1
    if(verbose) message(paste("      - (!) ENTREZ_ID values are missed for", 
                  length(df_rri$protein_id[which(is.na( df_rri$entrez_id[which(df_rri$is_contaminant == FALSE)] ))]), 
                  "protein ids"))
  }else{
    if(verbose) message("   + (+) All ENTREZ_ID ids available")
  }
  
  if( any(is.na(df_rri$redundant_ids)) ){
    # ic_rii <- ic_rii + 1
    if(verbose) message("      - (!) Some REDUNDANT_IDS values are missed (it should be fine)")
  }else{
    if(verbose) message("   + (+) All REDUNDANT_IDS values available")
  }
  
  if( any(is.na(df_rri$sequence)) ){
    ic_rii <- ic_rii + 1
    if(verbose) message("      - (-) ERROR: Some SEQUENCE values are missed")
  }else{
    if(verbose) message("   + (+) All SEQUENCE values available")
  }
  
  if("is_contaminant" %in% colnames(df_rri)){
    if( any(is.na(df_rri$is_contaminant)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) ERROR: Some IS_CONTAMINANT values are missed")
    }else{
      if(verbose) message("   + (+) All is_contaminant values available")
    }
    
    if(is.logical(df_rri$is_contaminant)){
      if(verbose) message("   + (+) IS_CONTAMINANT is a logical variable")
    }else{
      ic_rii <- ic_rii + 1
      if(verbose) message("   + (+) IS_CONTAMINANT is NOT a logical variable (all values should be TRUE/FALSE)")
    }
  }else{
    ic_rii <- ic_rii + 1
    if(verbose) message("      - (-) ERROR: IS_CONTAMINANT column is not available")
  }
  
  if("peptide_score" %in% colnames(df_rri)){
    
    if( any(is.na(df_rri$peptide_score)) ){
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) ERROR: Some PEPTIDE_SCORE values are missed")
    }else{
      if(verbose) message("   + (+) All PEPTIDE_SCORE values available")
    }
    
    if(is.numeric(df_rri$peptide_score)){
      if(verbose) message("   + (+) All PEPTIDE_SCORE values are numeric")
    }else{
      ic_rii <- ic_rii + 1
      if(verbose) message("      - (-) ERROR: PEPTIDE_SCORE is not numeric!")
    }
  }
  
  if( any(is.na(df_rri$organism_name)) ){
    ic_rii <- ic_rii + 1
    if(verbose) message("      - (-) ERROR: Some 	ORGANISM_NAME values are missed")
  }else{
    if(verbose) message("   + (+) All ORGANISM_NAME values available")
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
  tmt_channel <- NULL

  if("tmt11_channel" %in% colnames(df_vm)){
    required_columns <- c("vial_label", "tmt_plex", "tmt11_channel") 
    tmt_channel <- "tmt11_channel"
  }else if("tmt16_channel" %in% colnames(df_vm)){
    required_columns <- c("vial_label", "tmt_plex", "tmt16_channel") 
    tmt_channel <- "tmt16_channel"
  }else{
    if(verbose) message("      - (-) `tmt[11|16]_channel` column not found")
    ic_vm <- 3
    return(ic_vm)
  }
  
  if(!all( required_columns %in% colnames(df_vm))){
    if(verbose) message("      - (-) CRITICAL column(s) missed: ", appendLF = FALSE)
    if(verbose) message(paste(required_columns[!required_columns %in% colnames(df_vm)], collapse = ", "))
    ic_vm <- 3
    return(ic_vm)
  }
  
  if("tmt11_channel" %in% colnames(df_vm)){
    df_vm$vial_label <- gsub(" ", "", df_vm$vial_label)
    valid_channels <- c("126C", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C")
  }else if("tmt16_channel" %in% colnames(df_vm)){
    df_vm$vial_label <- gsub(" ", "", df_vm$vial_label)
    valid_channels <- c("126C", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N")
  }

  ic_vm <- 0

  plexes <- unique(df_vm$tmt_plex)

  for(p in plexes){
    temp_plex <- df_vm[which(df_vm$tmt_plex == p),]
    if(all(valid_channels %in% temp_plex[[tmt_channel]])){
      if(verbose) message("   + (+) All ", tmt_channel," channels are valid in plex ", paste(p))
    }else{
      if(tmt_channel == "tmt11_channel"){
        # Broad exception: they removed one plexed due to some issues. This is an exception to handle it
        valid_channels <- c("126C", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "131N", "131C")
        if(all(valid_channels %in% temp_plex$tmt11_channel)){
          if(verbose) message("   + ( ) Channel 130C missed in ", paste(p))
          ic_vm <- ic_vm + 1
        }else{
          if(verbose) message("      - (-) TMT Channels are missed in ", paste(p))
          ic_vm <- ic_vm + 1
        }
      }
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

  if( any(duplicated(df_vm$vial_label)) ){
    if(verbose) message("      - (-) VIAL_METADATA: duplications detected in the `vial_label` column: FAIL")
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
    # v_m$vial_label <- gsub(" ", "", v_m$vial_label)
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
    
    # Rename
    if( !is.null(all_samples) ){
      required_columns <- get_required_columns(isPTM = isPTM,
                                               prot_file = "rii")
      required_columns <- c(required_columns, sort(all_samples))

      if( all(required_columns %in% colnames(peprii)) ){
        peprii <- subset(peprii, select = required_columns)
      }else{
        what_is_missed <- required_columns[!(required_columns %in% colnames(peprii))]
        stop("RII: required columns from vial_label are not available in RII file. \nMissed columns: ", paste(what_is_missed, collapse = ","))
      }
    }

  }else{
    if(verbose) message("      - (-) {results_RII-peptide} file not available")
    ic <- ic + 10
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

    if( !is.null(all_vial_labels) ){
      required_columns <- get_required_columns(isPTM = isPTM,
                                               prot_file = "ratio")
      required_columns <- c(required_columns, all_vial_labels)

      if( all(required_columns %in% colnames(ratior)) ){
        ratior <- subset(ratior, select = required_columns)
      }else{
        stop("RATIO: required columns from vial_label are not available in RATIO file")
      }
    }
  }else{
    if(verbose) message("      - (-) {results_ratio.txt} file not available")
    ic <- ic + 10
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
#' @param dmaqc_shipping_info (char) DMAQC file with shipping information. If not provided
#' (default), then DMQAC validation is not peformed (only done at the BIC)
#' @param dmaqc_phase2validate (char) Provide phase to validate. This argument
#' is not required since it should be extracted from the input folder or from the 
#' new required file `metadata_phase.txt`. However, if this argument is provided,
#' it will take priority (and the phase from the input folder and the 
#' `metadata_phase.txt` will be ignored). Examples
#' - Folder with `PASS1A-06`: type either `PASS1A-06` or leave it `NULL`
#' - Both `PASS1A-06` and `PASS1C-06`: type `PASS1A-06|PASS1C-06`
#' - Only `PASS1C-06`: type `PASS1C-06`
#' @param f_proof (char) print out pdf with charts including:
#' - Reported Ion Intensity boxplot distribution and percentage of NA values per sample
#' - Ratio: ratio boxplot distribution and percentage of NA values per samples
#' @param out_qc_folder (char) if f_proof is TRUE, then a folder must be provided
#' @param return_n_issues (logical) if `TRUE` (default) returns the number of issues
#' @param full_report (logical) if `FALSE` (default) it returns only the
#' total number of issues. If `TRUE` returns the details of the number of issues (by
#' group of files, e.g., results, metadata_metabolites, etc)
#' @param printPDF (logical) if `TRUE` (default print plots to pdf)
#' @param check_only_results (logical) if `TRUE`, only validates results (default `FALSE`)
#' @param verbose (logical) `TRUE` (default) prints QC details.
#' @return (data.frame) Summary of issues
#' @export
validate_proteomics <- function(input_results_folder,
                                isPTM,
                                cas,
                                dmaqc_shipping_info = NULL,
                                dmaqc_phase2validate = NULL,
                                f_proof = FALSE,
                                out_qc_folder = NULL,
                                return_n_issues = TRUE,
                                full_report = FALSE,
                                printPDF = TRUE,
                                verbose = TRUE,
                                check_only_results = FALSE){
  
  percent_coverage <- NULL

  if(any(missing(input_results_folder) |
         missing(isPTM) |
         missing(cas)))
    stop("One (or many) of the required arguments missed.
        Please, check the help for this function to find out more")

  input_results_folder <- normalizePath(input_results_folder, winslash = "/")

  # Validate folder structure-----
  validate_cas(cas = cas)
  processfolder <- validate_processFolder(input_results_folder)
  assay <- validate_assay(input_results_folder)
  phase <- validate_phase(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)
  batch_folder <- validate_batch(input_results_folder)
  
  if(verbose) message("# PROTEOMICS QC report\n\n")
  
  # Set phase-----
  dmaqc_phase2validate <- set_phase(input_results_folder = input_results_folder, 
                                     dmaqc_phase2validate = dmaqc_phase2validate)
  
  phase2file <- gsub("\\|", "_", dmaqc_phase2validate)

  # Print out proofs----
  if(f_proof){
    if(!is.null(out_qc_folder)){
      if(!dir.exists(file.path(out_qc_folder))){
        dir.create(file.path(out_qc_folder), recursive = TRUE)
      }
    }else{
      out_qc_folder <- getwd()
    }

    output_prefix <- paste0(cas, ".", tolower(phase2file), ".", tissue_code, ".",tolower(assay), ".", tolower(processfolder))
  }

  input_folder_short <- regmatches(input_results_folder, regexpr("(PASS|HUMAN).*RESULTS_[0-9]{8}", input_results_folder))
  if( purrr::is_empty(input_folder_short) ){
    if(verbose) message("\nThe RESULTS_YYYYMMDD folder full path is not correct. Example:")
    if(verbose) message("/full/path/to/folder/PASS1A-06/T66/RPNEG/BATCH1_20190822/RESULTS_20200308")
    stop("Input folder not according to guidelines")
  }

  if(verbose) message("+ Site: ", toupper(cas))
  if(verbose) message("+ Folder: `", paste0(input_folder_short),"`")

  # COUNTS----
  ic <- 0 # Critical issues
  ic_rii <- NA # RII issues
  ic_rr <- NA # ratio results
  ic_vm <- NA # vial metadata
  ic_man <- 0 # manifest

  if(missing(dmaqc_shipping_info)){
    dmaqc_shipping_info = NULL
  }

  if(is.null(dmaqc_shipping_info)){
    ic_vl <- "missed"
  } else{
    ic_vl <- NA
  }

  # VIAL METADATA-----

  if(verbose) message("\n## VIAL METADATA\n")

  lista <- open_file(input_results_folder = input_results_folder,
                     filepattern = "vial_metadata.txt",
                     verbose = verbose)

  f_vm <- lista$flag

  all_vial_labels <- NULL

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
    
    # # PlexedPiper temporal error
    # if( all(c("Ref_S1", "Ref_S2", "Ref_S3", "Ref_S4", "Ref_S5", "Ref_S6") %in% colnames(peprii)) ){
    #   Ref_S1=Ref_S2=Ref_S3=Ref_S4=Ref_S5=Ref_S6=NULL
    #   peprii <- rename(peprii, Ref_A=Ref_S1)
    #   peprii <- rename(peprii, Ref_B=Ref_S2)
    #   peprii <- rename(peprii, Ref_C=Ref_S3)
    #   peprii <- rename(peprii, Ref_D=Ref_S4)
    #   peprii <- rename(peprii, Ref_E=Ref_S5)
    #   peprii <- rename(peprii, Ref_F=Ref_S6)
    # }
    
    ic_rii <- check_rii_proteomics(df_rri = peprii,
                                   isPTM = isPTM,
                                   f_proof = f_proof,
                                   output_prefix = output_prefix,
                                   return_n_issues = TRUE,
                                   out_qc_folder = out_qc_folder,
                                   printPDF = printPDF,
                                   verbose = verbose)

    if( is.null(ic_rii) ){
      f_rii <- FALSE
      ic_rii <- NA
      ic <- ic + 1
    }else{
      if(f_proof & f_vm){
        # print proof function here
        # Only if vial_label metadata is available
        proteomics_plots_rii(all_vial_labels = all_vial_labels,
                             all_samples = all_samples,
                             peprii = peprii,
                             isPTM = isPTM,
                             v_m = v_m,
                             out_qc_folder = out_qc_folder,
                             output_prefix = output_prefix,
                             printPDF = printPDF,
                             verbose = verbose)
      } # print plots
    }
  }else{
    if(verbose) message("      - (-) {results_RII-peptide} file not available")
    ic <- ic + 10
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

        if(verbose) message("   + (+) PLOTS RATIO------------------")
        
        # Get the total number of samples to customize the plots. If larger than 200, 
        # prepare for large plots
        sn <- length(unique(all_samples))
        # Set a limit for which remove labels 
        sn_limit <- 200

        if( !is.null(all_vial_labels) ){
          required_columns <- get_required_columns(isPTM = isPTM,
                                                   prot_file = "ratio")
          required_columns <- c(required_columns, all_vial_labels)

          if( all(required_columns %in% colnames(ratior)) ){
            if(isPTM){
              r_c <- c("ptm_id", all_vial_labels)
              fratior <- subset(ratior, select = r_c)
              ratior_long <- fratior %>% tidyr::pivot_longer(cols = -c(ptm_id),
                                                            names_to = "vial_label",
                                                            values_to = "ratio_values")
            }else{
              r_c <- c("protein_id", all_vial_labels)
              fratior <- subset(ratior, select = r_c)
              ratior_long <- fratior %>% tidyr::pivot_longer(cols = -c(protein_id),
                                                            names_to = "vial_label",
                                                            values_to = "ratio_values")
            }

            if(f_vm){
              
              # Plotting ratio distributions-----
              if(verbose) message("       - (p) Plotting ratio distributions")

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

              # Plotting NA values-----
              if(verbose) message("       - (p) Plotting NA percentage in ratio results")
              p_na_ratior <- ratior[required_columns] %>%
                inspectdf::inspect_na() %>%
                dplyr::arrange(match(col_name, colnames(ratior))) %>%
                inspectdf::show_plot() +
                ylim(0, 100) + theme_linedraw() +
                theme(axis.text.x = element_text(angle = 90,
                                                 hjust = 1,
                                                 vjust = 0.5,
                                                 size = 8))
              
              # Plot unique IDs-----
              if(verbose) message("       - (p) Plotting unique ids in ratio")
              if(isPTM){
                key_id <- "ptm_id"
              }else{
                key_id <- "protein_id"
              }
              
              # uid <- ratior_long %>% 
              #   group_by(across(all_of(c(key_id, "vial_label", "tmt_plex")))) %>% 
              #   summarise(total_rii = sum(ratio_values, na.rm = FALSE), .groups = 'drop')
              uid <- ratior_long %>% 
                group_by(across(all_of(c(key_id, "vial_label", "tmt_plex")))) %>% 
                summarise(total_rii = ratio_values, .groups = 'drop')
              
              uid2 <- uid[which(!is.na(uid$total_rii)),]
              uid3 <- unique(uid2[c(key_id, "vial_label", "tmt_plex")]) %>% 
                count(vial_label, tmt_plex)
              
              puid1 <- ggplot(uid3, aes(x = reorder(vial_label, n), y = n, fill = tmt_plex)) + 
                geom_bar(stat = "identity") + 
                theme_linedraw() +
                theme(
                  axis.text.x = element_text(
                    angle = 90,
                    hjust = 1,
                    vjust = 0.5,
                    size = 8
                  ),
                  legend.position = "none"
                ) +
                geom_text(
                  aes(label = n),
                  # vjust = -0.5,
                  hjust = 1,
                  size = 2.7,
                  angle = 90
                ) +
                ggtitle("Ratio: Unique IDs in samples") +
                facet_wrap(~ tmt_plex, scales = "free") +
                xlab("Vial Labels")
              
              puid2 <- ggplot(uid3, aes(x = reorder(vial_label, n), y = n, fill = tmt_plex)) + 
                geom_bar(stat = "identity",
                         na.rm = TRUE) +
                theme_linedraw() +
                theme(
                  axis.text.x = element_text(
                    angle = 90,
                    hjust = 1,
                    vjust = 0.5,
                    size = 8
                  )
                ) +
                geom_text(
                  aes(label = n),
                  # vjust = -0.5,
                  hjust = 1,
                  size = 2.7,
                  angle = 90
                ) +
                ggtitle("Ratio: Unique IDs in samples") +
                xlab("Vial Labels")
              
              # Plotting protein coverage-----
              if(!isPTM){
                if(verbose) message("       - (p) Plotting protein coverage")
                ppp <- ggplot(ratior, aes(x=percent_coverage)) +
                  geom_histogram( binwidth=2, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
                  ggtitle("Overall Protein Identification Coverage") +
                  theme_light() +
                  theme(plot.title = element_text(size=14)) + 
                  labs(x = "Percent Coverage", y = "Number of Proteins")
                
                pbp <- ggplot(ratior, aes(y = percent_coverage)) + 
                  geom_boxplot(col="#69b3a2") +
                  ggtitle("") +
                  theme_light() +
                  theme(plot.title = element_text(size=14)) + 
                  theme(axis.title.x=element_blank())
              }
              
              # + scale_fill_brewer(palette="Reds")
              
              if(is.null(out_qc_folder)){
                out_plot_ratdist <- paste0(output_prefix,"-qc-ratio-distribution.pdf")
              }else{
                out_plot_ratdist <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-ratio-distribution.pdf"))
              }

              if(printPDF){
                if(sn > 800){
                  pdf(out_plot_ratdist, width = 40, height = 8)
                }else if(sn <= 800 & sn > 200){
                  pdf(out_plot_ratdist, width = 22, height = 8)
                }else{
                  pdf(out_plot_ratdist, width = 14, height = 8)
                }
              }
                print(pisr)
                print(puid1)
                print(puid2)
                print(p_na_ratior)
                if(!isPTM) grid.arrange(arrangeGrob(ppp, pbp, ncol = 2))
              if(printPDF) garbage <- dev.off()
            }
          }else{
            if(verbose) message("      - (-) Required columns are missed in ratio file (based on vial_label). Print proof is not possible")
            ic_rr <- ic_rr + 1
          }
        }else{
          if(verbose) message("      - (-) Vial labels are NULL: not possible to generate proof")
          ic_rr <- ic_rr + 1
        }
      }# if proof
    } #print plots
  }else{
    if(verbose) message("      - (-) {results_ratio.txt} file not available")
    ic <- ic + 10
  }

  # MANIFEST----
  if(check_only_results){
    f_man <- FALSE
    ic_man <- 0
  }else{
    
    if(verbose) message("\n## MANIFEST\n")
    
    file_manifest <- list.files(normalizePath(batch_folder),
                                pattern="file_manifest",
                                ignore.case = TRUE,
                                full.names=TRUE,
                                recursive = TRUE)
    
    if(length(file_manifest) == 0){
      f_man <- FALSE
      ic_man <- ic_man + 1
    }else if(length(file_manifest) >= 1){
      file_manifest <- file_manifest[length(file_manifest)]
      f_man <- TRUE
    }
    
    if(f_man){
      
      # load all files in manifest:
      f_list = list()
      for (f in file_manifest ){
        f_list[[f]] = read.csv(f)
        f_list[[f]]$file <- f
      }
      
      manifest <- bind_rows(f_list)
      
      manifest$file_name <- basename(manifest$file_name)
      
      mani_columns <- c("file_name", "md5")
      if( all(mani_columns %in% colnames(manifest)) ){
        if(verbose) message("   + (+) <file_name, md5> columns available in manifest file")
        if(f_rr){
          ratio_file <- manifest$file_name[grepl("ratio.txt", manifest$file_name)]
          if( length(ratio_file) == 1 ){
            if( file.exists(file.path(input_results_folder, ratio_file)) ){
              if(verbose) message("   + (+) RATIO file included")
            }else{
              if(verbose) message("      - (-) The ratio file name provided in manifest is not found in folder")
              ic_man <- ic_man + 1
            }
          }else if (length(ratio_file) > 1){
            if(verbose) message("      - (-) More than one ratio.txt files available in the folder")
            ic_man <- ic_man + 1
          }else{
            # This should never happen
            if(verbose) message("      - (-) No ratio.txt file available in the manifest file")
            ic_man <- ic_man + 1
          }
        }
        
        if(f_rii){
          rii_file <- manifest$file_name[grepl("RII-peptide.txt", manifest$file_name, ignore.case = TRUE)]
          if(length(rii_file) == 1){
            if(file.exists(file.path(input_results_folder, rii_file))){
              if(verbose) message("   + (+) RII file included")
            }else{
              if(verbose) message("      - (-) RII file name provided in manifest is not found in folder")
              ic_man <- ic_man + 1
            }
          }else if (length(rii_file) > 1){
            if(verbose) message("      - (-) More than one RII-peptide.txt files available in the folder")
            ic_man <- ic_man + 1
          }else{
            # This should never happen
            if(verbose) message("      - (-) No RII-peptide.txt available in the manifest file")
            ic_man <- ic_man + 1
          }
          
        }
        
        if(f_vm){
          vm_file <- manifest$file_name[grepl("vial_metadata", manifest$file_name, ignore.case = TRUE)]
          if( length(vm_file) == 1 ){
            if(file.exists(file.path(input_results_folder, vm_file))){
              if(verbose) message("   + (+) VIAL_METADATA file included")
            }else{
              if(verbose) message("      - (-) vial_metadata name provided in manifest is not found in folder")
              ic_man <- ic_man + 1
            }
          }else if( length(vm_file) > 1 ) {
            if(verbose) message("      - (-) VIAL_METADATA file is not included in manifest file")
            ic_man <- ic_man + 1
          }else{
            if(verbose) message("      - (-) VIAL_METADATA file is not included in manifest file")
            ic_man <- ic_man + 1
          }
        }
        
        if( any(is.na(manifest$md5)) ){
          if(verbose) message("      - (-) MD5 column contains NA values")
          ic_man <- ic_man + 1
        }
        
      }else{
        if(verbose) message("      - (-) MANIFEST ERROR: Missing columns (must contain <file_name> and <md5>)")
        ic_man <- ic_man + 1
        ic <- ic + 1
      }
    }else{
      message("      - (-) ERROR: manifest file not found")
      ic_man = 6
    }
  } #check_only_results?
  
  if(ic_man > 0){
    ic <- ic + ic_man
  }

  # INTER-FILE VALIDATION----

  if(verbose) message("\n## INTER-file validations\n")

  if(f_vm & f_rii & f_rr){

    if( all(all_vial_labels %in% colnames(peprii)) ){
      if(verbose) message("   + (+) All vial labels samples available in RII file")
    }else{
      ic <- ic + 1
      if(verbose) message("      - (-) CRITICAL: Some vial labels are not available in RII file: FAIL")
    }

    if( all(all_vial_labels %in% colnames(ratior)) ){
      if(verbose) message("   + (+) All vial labels available in RATIO results file")
    }else{
      ic <- ic + 1
      if(verbose) message("      - (-) CRITICAL: Some vial labels not available in RATIO results file: FAIL")
    }

  }else{
    if(verbose) message("      - (-) CRITICAL: One or many of the files are not available: FAIL")
    ic <- ic + 3
  }

  # CHECK DMAQC----

  # Validate vial labels from DMAQC

  if(verbose) message("\n\n## DMAQC validation\n")
  failed_samples <- check_failedsamples(input_results_folder = input_results_folder, 
                                        verbose = verbose)

  if( is.na(ic_vl) ){
    if(f_vm){
      outfile_missed_viallabels <- paste0(cas, ".", tolower(phase2file), ".", tissue_code, ".",tolower(assay), ".", tolower(processfolder))
      
      ic_vl <- check_viallabel_dmaqc(vl_submitted = all_vial_labels,
                                     dmaqc_shipping_info = dmaqc_shipping_info,
                                     tissue_code = tissue_code,
                                     cas = cas,
                                     phase = dmaqc_phase2validate, 
                                     failed_samples = failed_samples,
                                     out_qc_folder = out_qc_folder,
                                     outfile_missed_viallabels = outfile_missed_viallabels,
                                     return_n_issues = TRUE,
                                     verbose = verbose)
    }else{
      ic <- ic + 1
      if(verbose) message("      - (-) DMAQC validation cannot be performed (no vial label data available)")
    }
  }


  # PRINT OUT RESULTS-----
  batchversion <- stringr::str_extract(string = input_results_folder, pattern = "BATCH.*_[0-9]+/RESULTS_[0-9]+")

  qc_date <- Sys.time()
  qc_date <- gsub("-", "", qc_date)
  qc_date <- gsub(" ", "_", qc_date)
  qc_date <- gsub(":", "", qc_date)
  t_name <- bic_animal_tissue_code$bic_tissue_name[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)]

  if(return_n_issues){
    total_issues <- sum(ic, ic_vm, ic_rr, ic_rii, ic_man, na.rm = TRUE)
    if(verbose) message("\nTOTAL NUMBER OF ISSUES: ", total_issues,"\n")
    if(full_report){
      reports <- data.frame(cas = cas,
                            phase= dmaqc_phase2validate,
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
#' @param folder_name (char) output folder name.
#' @param folder_root (char) absolute path to write the output folder. Default: current directory
#' @param version_file (char) file version number (v#.#)
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
                                      folder_name = "motrpac_release",
                                      folder_root = NULL,
                                      version_file = "v1.0",
                                      verbose = TRUE){

  # Get names from input_results_folder------
  assay <- validate_assay(input_results_folder)
  phase <- validate_phase(input_results_folder)
  tissue_code <- validate_tissue(input_results_folder)
  folder_phase <- tolower(phase)
  folder_tissue <- bic_animal_tissue_code$tissue_name_release[which(bic_animal_tissue_code$bic_tissue_code == tissue_code)]
  
  phase_metadata <- set_phase(input_results_folder = input_results_folder, 
                              dmaqc_phase2validate = NULL, 
                              verbose = FALSE)
  phase_details <- generate_phase_details(phase_metadata)

  assay_codes$cas_code <- NULL
  assay_codes <- unique(assay_codes)
  if(length(assay_codes$assay_code[which(assay_codes$submission_code == assay)]) == 1){
    folder_assay <- assay_codes$assay_code[which(assay_codes$submission_code == assay)]
  }else{
    stop("ASSAY code ", assay, " not available in < assay_codes >")
  }
  

  if(verbose) message("+ Writing out ", phase, " (phase-details: ", phase_details, ") ", tissue_code, " ", assay, " files", appendLF = FALSE)

  if( grepl("PH", assay) | grepl("AC", assay) | grepl("UB", assay)){
    isPTM = TRUE
  }else{
    isPTM = FALSE
  }
  
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

  output_folder <- file.path(folder_root, folder_name, folder_phase, "proteomics-untargeted", folder_tissue, folder_assay)

  if(!dir.exists(file.path(output_folder))){
    dir.create(file.path(output_folder), recursive = TRUE)
  }

  file_name_shared <- paste0("motrpac_",
                             phase_details, "_",
                             folder_tissue, "_",
                             folder_assay)


  # Create and write FILES-----
  prot_rii <- file.path(output_folder, paste0(file_name_shared,"_rii-results_", version_file, ".txt"))
  prot_ratio <- file.path(output_folder, paste0(file_name_shared,"_ratio-results_", version_file, ".txt"))
  vial_metadata <- file.path(output_folder, paste0(file_name_shared,"_vial-metadata_", version_file, ".txt"))

  write.table(prot_dfs$peprii, prot_rii, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(prot_dfs$ratior, prot_ratio, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(prot_dfs$v_m, vial_metadata, row.names = FALSE, sep = "\t", quote = FALSE)

  if(verbose) message("...done!")
}

# @title Get proteomics required columns
get_required_columns <- function(isPTM, prot_file){
  if(prot_file == "ratio"){
    both_required <- c("protein_id", 
                       "gene_symbol", 
                       "entrez_id", 
                       "redundant_ids", 
                       "organism_name", 
                       "is_contaminant")
    if(isPTM){
      required_columns <- c(both_required, 
                            "ptm_id", 
                            "confident_score", 
                            "confident_site", 
                            "flanking_sequence", 
                            "ptm_score")
      
    }else{
      required_columns <- c(both_required, 
                            "num_peptides", 
                            "percent_coverage", 
                            "protein_score")
    }
  }else if(prot_file == "rii"){
    both_required <- c("protein_id", 
                       "redundant_ids",
                       "is_contaminant",
                       "peptide_score", 
                       "sequence", 
                       "gene_symbol", 
                       "entrez_id", 
                       "organism_name")
    if(isPTM){
      required_columns <- c(both_required,
                            "ptm_id", 
                            "ptm_peptide", 
                            "confident_score", 
                            "confident_site")
    }else{
      required_columns <- c(both_required)
    }
  }
  return(required_columns)
}

