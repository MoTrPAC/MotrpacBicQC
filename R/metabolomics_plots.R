
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot Basic Metabolomics QC charts
#'
#' @description Plot intensity distributions, number of unique ids, NA values 
#' per sample
#' @param results (df) metabolomics results (both named and unnamed already merged, if untargeted)
#' @param m_s_n (df) metadata samples named
#' @param out_qc_folder (char) output qc folder (it creates the folder if it doesn't exist)
#' @param output_prefix (char) prefix for the file name output (pdf file)
#' @param printPDF (logical) `TRUE` (default) prints pdf file
#' @param verbose (logical) `TRUE` (default) shows messages
#' @export
plot_basic_metabolomics_qc <- function(results, 
                                       m_s_n,
                                       out_qc_folder = NULL,
                                       output_prefix,
                                       printPDF = TRUE,
                                       verbose = TRUE){

  metabolite_name = id_type = sample_id = sample_order = intensity = sample_type = NULL
  
  if(verbose) message("   + (+) QC PLOTS ------------------")
  
  results_long <- results %>% tidyr::pivot_longer(cols = -c(metabolite_name, id_type),
                                                  names_to = "sample_id",
                                                  values_to = "intensity")
  
  if(verbose) message("       - (p) Plot intensity distributions")
  
  results_long <- merge(m_s_n, results_long, by = c("sample_id"))
  
  results_long$sample_id <- as.character(results_long$sample_id)
  results_long$sample_id <- as.factor(results_long$sample_id)
  results_long$sample_type <- as.factor(results_long$sample_type)
  
  piseso <- ggplot2::ggplot(results_long,
                            aes(
                              # x = reorder(sample_id, log2(intensity), FUN = median, na.rm = TRUE),
                              x = reorder(sample_id, sample_order),
                              y = log2(intensity),
                              fill = sample_type)) +
    geom_boxplot(na.rm = TRUE) +
    theme_linedraw() +
    theme(
      # axis.text.x = element_text(angle = 90,
      # hjust = 1,
      # vjust = 0.5,
      # size = 8),
      axis.title.x = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 90),
      legend.position="top") + #legend.position = "none"
    labs(title = "Intensity distribution / Unique IDs",
         subtitle = paste(output_prefix), 
         # x = "sample_id",
         # caption = "Sort by injection order",
         x = "",
         y = "log2(intensity)"
    )
  

  piseio <- ggplot2::ggplot(results_long,
                            aes(x = reorder(sample_id, sample_order),
                              # x = reorder(sample_id, log2(intensity), FUN = median, na.rm = TRUE),
                                y = log2(intensity),
                                fill = sample_type)) +
    geom_boxplot(na.rm = TRUE) +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 8)) + #legend.position = "none"
    labs(x = "sample_id", y = "log2(intensity)") +
    labs(title = "Intensity distribution",
         subtitle = paste(output_prefix),
         caption = "Sort by injection order")
  
  # Get data ready to count IDs
  if(verbose) message("       - (p) Plot ID counts")
  uid <- results_long %>%
    group_by(across(all_of(c("metabolite_name", "sample_id", "sample_type", "sample_order", "id_type")))) %>%
    summarise(total_intensity = intensity, .groups = 'drop')
  uid2 <- uid[which(!is.na(uid$total_intensity)),]
  uid3 <- unique(uid2[c("metabolite_name", "sample_id", "sample_type", "sample_order", "id_type")]) %>%
    count(sample_id, sample_type, sample_order, id_type)
  
  puid1 <- ggplot(uid3, aes(x = reorder(sample_id, sample_order), y = n, fill = sample_type)) +
    geom_bar(stat = "identity",
             position="stack",
             na.rm = TRUE) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      # legend.position = "none",
      axis.text.y = element_text(angle = 90)
    ) +
    geom_text(
      aes(label = n),
      # vjust = -0.5,
      hjust = 1,
      size = 2.7,
      angle = 90
    ) +
    labs(title = "Unique IDs in samples", 
         subtitle = paste(output_prefix),
         caption = "Scales are different",
         x = "Sample ID") +
    facet_grid(rows = vars(id_type), scales = "free")
  
  uid4 <- unique(uid2[c("metabolite_name", "sample_id", "sample_type", "sample_order")]) %>%
    count(sample_id, sample_type, sample_order)
  
  puid2 <- ggplot(uid4, aes(x = reorder(sample_id, sample_order), y = n, fill = sample_type)) +
    geom_bar(stat = "identity",
             na.rm = TRUE) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      axis.text.y = element_text(angle = 90),
      legend.position = "none"
    ) +
    geom_text(
      aes(label = n),
      # vjust = -0.5,
      hjust = 1,
      size = 2.7,
      angle = 90
    ) +
    # ggtitle("Unique IDs in samples") +
    labs(caption = "Sort by injection order") +
    xlab("Sample IDs")
  
  # Plot NA values------
  if(verbose) message("       - (p) Plot NA values")
  
  p_na_peprii <- results %>%
    inspectdf::inspect_na() %>%
    dplyr::arrange(match(col_name, colnames(results))) %>%
    inspectdf::show_plot() +
    ylim(0, 100) + theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 8)) +
    labs(title = "Prevalence of NAs", 
         subtitle = paste(output_prefix))
  
  
  if(!is.null(out_qc_folder)){
    if(!dir.exists(file.path(out_qc_folder))){
      dir.create(file.path(out_qc_folder), recursive = TRUE)
    }
  }else{
    out_qc_folder <- getwd()
  }

  out_plot_dist <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-basic-plots.pdf"))
  
  if(printPDF) pdf(out_plot_dist, width = 12, height = 8)
  print(puid1)
  print(piseio)
  gridExtra::grid.arrange(piseso, puid2, ncol = 1, heights = c(2, 1))
  print(p_na_peprii)
  if(printPDF) garbage <- dev.off()
}


