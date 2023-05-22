
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title QC Plot of proteomics reporter ion intensity data
#'
#' @description check whether results file is following guidelines
#' @param all_samples (vector) all sample ids
#' @param all_vial_labels (vector) only vial_labels
#' @param peprii (char) Reporter Ion Intensity df
#' @param isPTM (logical) Is a PTM assay
#' @param v_m (df) sample metadata
#' @param out_qc_folder (char) if `f_proof is TRUE`, a folder path can be provided
#' (otherwise print in current working directory)
#' @param output_prefix (char) provide a prefix for the output name
#' @param printPDF (logical) if `TRUE` (default print plots to pdf)
#' @param verbose (logical) `TRUE` (default) shows messages
#' @return (int) number of issues identified
#' @examples {
#' check_results(r_m = results_named, m_s = metadata_sample_named, m_m = metadata_metabolites_named)
#' }
#' @export
proteomics_plots_rii <- function(all_vial_labels,
                                 all_samples, 
                                 peprii,
                                 isPTM,
                                 v_m,
                                 out_qc_folder = NULL,
                                 output_prefix,
                                 printPDF,
                                 verbose){
                             
  if(verbose) message("   + (+) PLOTS RII------------------")
  
  # Get the total number of samples to customize the plots. If larger than 200, 
  # prepare for large plots
  sn <- length(unique(all_samples))
  # Set a limit for which remove labels 
  sn_limit <- 200
  
  
  if( !is.null(all_vial_labels) ){
    required_columns <- get_required_columns(isPTM = isPTM,
                                             prot_file = "rii")
    required_columns <- c(required_columns, all_samples)
    
    if( all(required_columns %in% colnames(peprii)) ){
      # Check distributions
      
      if(isPTM){
        r_c <- c("ptm_peptide", all_samples)
        fpeprii <- subset(peprii, select = r_c)
        peptides_long <- fpeprii %>% tidyr::pivot_longer(cols = -c(ptm_peptide),
                                                         names_to = "vial_label",
                                                         values_to = "ri_intensity")
      }else{
        r_c <- c("protein_id", "sequence", all_samples)
        fpeprii <- subset(peprii, select = r_c)
        peptides_long <- fpeprii %>% tidyr::pivot_longer(cols = -c(protein_id, sequence),
                                                         names_to = "vial_label",
                                                         values_to = "ri_intensity")
      }
      
      # Plot Intensity distribution-----
      if(verbose) message("       - (p) Plot intensity distributions")
      
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
      
      
      # Plot NA values------
      if(verbose) message("       - (p) Plot NA values")
      
      p_na_peprii <- peprii[required_columns] %>%
        inspectdf::inspect_na() %>%
        dplyr::arrange(match(col_name, colnames(peprii))) %>%
        inspectdf::show_plot() +
        ylim(0, 100) + theme_linedraw() +
        theme(axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5,
                                         size = 8))
      
      # Plot unique IDs-----
      if(verbose) message("       - (p) Plot Unique IDs")
      
      if(isPTM){
        key_id <- "ptm_peptide"
      }else{
        key_id <- "protein_id"
      }
      
      uid <- peptides_long %>% 
        group_by(across(all_of(c(key_id, "vial_label", "tmt_plex")))) %>% 
        dplyr::reframe(total_rii = ri_intensity)
      
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
        ggtitle("RII: Unique IDs in samples") +
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
        ggtitle("RII: Unique IDs in samples") + 
        xlab("Vial Labels")
      
      if(is.null(out_qc_folder)){
        out_plot_dist <- paste0(output_prefix,"-qc-rii-distribution.pdf")
      }else{
        out_plot_dist <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-rii-distribution.pdf"))
      }
      
      if(printPDF){
        if(sn > 800){
          pdf(out_plot_dist, width = 40, height = 8)
        }else if(sn <= 800 & sn > 200){
          pdf(out_plot_dist, width = 22, height = 8)
        }else{
          pdf(out_plot_dist, width = 14, height = 8)
        }
      }
      print(pise)
      print(puid1)
      print(puid2)
      print(p_na_peprii)
      if(printPDF) garbage <- dev.off()
    }else{
      if(verbose) message("      - (-) Not all required columns available in RII: print proof is not possible")
      required_columns[!(required_columns %in% colnames(peprii))]
      ic <- ic + 1
    }
  } # !is.null(all_vial_labels)
}
