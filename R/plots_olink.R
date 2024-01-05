
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot Basic OLINK QC charts
#'
#' @description Plot value distributions, number of unique ids, NA values
#' @param results (df) olink results merged with proteins metadata
#' @param results_long (df) olink results in long format and merged with 
#' sample metadata
#' @param out_qc_folder (char) output qc folder (it creates the folder if it doesn't exist)
#' @param output_prefix (char) prefix for the file name output (pdf file)
#' @param printPDF (logical) `TRUE` (default) prints pdf file
#' @param verbose (logical) `TRUE` (default) shows messages
#' @export
plot_basic_olink_qc <- function(results, 
                                results_long,
                                out_qc_folder = NULL,
                                output_prefix,
                                printPDF = TRUE,
                                verbose = TRUE){
                               

  olink_id = panel_name = sample_id = sample_id_ordered = value = sample_type = olink_count = plate_id = NULL
  
  if(verbose) message("   + (+) QC PLOTS ------------------")
  
  if(verbose) message("       - (p) Plot value distributions")
  
  # Get the total number of samples to customize the plots. If larger than 200, 
  # prepare for large plots
  
  sn <- length(unique(results_long$sample_id))
  # Set a limit for which remove labels 
  sn_limit <- 200
  
  # Density plots------ 
  
  if(verbose) message("       - (p) Density distributions")
  mu <- results_long %>%
    group_by(sample_id) %>%
    dplyr::reframe(grp.mean = mean(value))
  
  den1 <- ggplot(data = results_long, aes(x = value, color = sample_type)) + 
    geom_density(na.rm = TRUE) + 
    theme_light() +
    labs(title = "Density distribution by sample type", 
         subtitle = paste(output_prefix)) 
  
  den2 <- ggplot(data = results_long, aes(x = value, 
                                          color = sample_id)) +
    geom_density(na.rm = TRUE) + 
    theme_light() +
    theme(legend.position="none") +
    labs(title = "Density distribution by sample id", 
         subtitle = paste(output_prefix))
  
  den3 <- ggplot(data = results_long, aes(x = value, 
                                          color = sample_id)) +
    geom_density(na.rm = TRUE) + 
    theme_minimal() +
    theme(legend.position="none") + 
    facet_wrap(~sample_type) +
    labs(title = "Density distributions", 
         subtitle = paste(output_prefix)) 
  
  # Boxplot distribution, as part of a grid----
  piseso <- ggplot(results_long, aes(x = sample_id_ordered, 
                                     y = value, 
                                     fill = plate_id)) +
    geom_boxplot(na.rm = TRUE) +
    theme_classic() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 90),
          legend.position = "top") +
    labs(title = "value distribution / Unique IDs", 
         subtitle = paste(output_prefix), 
         y = "value")
  
  # Boxplots: standalone----
  piseio <- 
    ggplot2::ggplot(results_long, aes(x = sample_id_ordered, 
                                      y = value, 
                                      fill = sample_type)) +
    geom_boxplot(na.rm = TRUE) +
    theme_classic() +
    labs(
      x = "sample_id", 
      subtitle = paste(output_prefix),
      caption = "Sort by injection order"
    ) + 
    theme(legend.position="top")

  if(sn > sn_limit){
    piseio <- piseio + theme(
      axis.text.x = element_blank())
  }else{
    piseio <- piseio + theme(
      axis.text.x = element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5,
                                   size = 6))
  }

  piseio <- piseio + labs(title = "Sample value distribution",
                          y = "value") 
  
  # Plot number counts
  if(verbose) message("       - (p) Plot ID counts")
  count_data <- results_long %>%
    group_by(plate_id, sample_id) %>%
    summarise(olink_count = n_distinct(olink_id), .groups = 'drop')
  
  count_data <- count_data %>%
    arrange(plate_id) %>%
    mutate(sample_id_ordered = factor(sample_id, levels = unique(sample_id)))
  
  puid1 <- ggplot(count_data, aes(x = sample_id_ordered, y = olink_count, fill = plate_id)) +
    geom_bar(stat = "identity", 
             position = "stack") +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    labs(title = "Counts of unique olink_id by Sample and Plate",
         x = "Samples",
         y = "Count of olink_id") +
    theme(axis.text.x = element_blank(), legend.position = "top")
  if(sn > sn_limit){
    puid1 <- puid1 + theme(axis.text.x = element_blank())
  }else{
    puid1 <- puid1 + theme(
      axis.text.x = element_text(angle = 90,
                                 hjust=0.95,vjust=0.2,
                                 size = 6)) +
      geom_text(
        aes(label = n),
        hjust = 1,
        size = 2,
        angle = 90
      ) 
  }
  
  # Plot number of ids per sample-----
  
  puid2 <-  ggplot(count_data, aes(x = sample_id_ordered, y = olink_count, fill = plate_id)) +
    geom_bar(stat = "identity", 
             position = "stack") +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    labs(x = "Samples",
         y = "Count of olink_id") +
    theme(axis.text.x = element_blank(), legend.position = "none")
  if(sn > sn_limit){
    puid2 <- puid2 + theme(axis.text.x = element_blank())
  }else{
    puid2 <- puid2 + theme(
      axis.text.x = element_text(angle = 90,
                                 hjust=0.95,vjust=0.2,
                                 size = 6)) +
      geom_text(
        aes(label = n),
        hjust = 1,
        size = 2,
        angle = 90
      ) + ylab(NULL)
  }
  
  # Plot proportion of ids------
  if(verbose) message("       - (p) Plot Sample counts and ID numbers/proportions")
  
  # Plot Total number of samples per sample_type
  tns <- unique(results_long[c("sample_id", "sample_type")]) %>% group_by(sample_type) %>% count(sample_type)
  ptns <- ggplot(tns, aes(x = sample_type, y = n, fill = sample_type)) +
    geom_bar(stat = "identity") +
    labs(title = "Number of samples by sample type", 
         subtitle = paste(output_prefix), 
         y = "", x = "") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14),
          ) + 
    geom_text(
      aes(label = n),
      hjust = 0.5,
      vjust = -0.2,
      size = 4,
      angle = 0
    ) +
    scale_fill_brewer(palette = "Dark2")
  
  
  ppids <- ggplot(results, aes(x = as.factor(panel_name), fill = as.factor(panel_name))) +
    geom_bar(aes(y = after_stat(count)/sum(after_stat(count)))) +
    geom_text(aes(y = after_stat(count)/sum(after_stat(count)), 
                  label = scales::percent(after_stat(count)/sum(after_stat(count)))), 
              stat = "count", vjust = -0.25) +
    scale_y_continuous(labels = percent) +
    labs(title = "Proportion of Features Identified (named/unnamed)", 
         subtitle = paste(output_prefix), 
         y = "", x = "") +
    theme_light() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18)) + 
    scale_fill_brewer(palette = "Dark2")
  
  pnids <- ggplot(results, aes(x = panel_name, fill = panel_name)) +
    geom_bar() +
    geom_text(stat = 'count',aes(label = after_stat(count), vjust = -0.2)) +
    labs(title = "Total Number of Features Identified (named/unnamed)", 
         subtitle = paste(output_prefix), y = "", x = "") +
    theme_light() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18)) + 
    scale_fill_brewer(palette = "Dark2")

  
  # Plot NA values------
  if(verbose) message("       - (p) Plot NA values")
  
  suppressWarnings(
    p_na_peprii <- results %>%
      inspectdf::inspect_na() %>% 
      dplyr::arrange(match(col_name, colnames(results))) %>% 
      inspectdf::show_plot(text_labels = FALSE) + ylim(0, 100) + 
      theme_classic() +
      labs(title = "Prevalence of NAs",
           subtitle = paste(output_prefix),
           y = "% of NAs in each sample"))
    
    if(sn > sn_limit){
      p_na_peprii <- p_na_peprii + 
        theme(axis.text.x = element_blank())
    }else{
      p_na_peprii <- p_na_peprii + 
        theme(axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5,
                                         size = 8))
    }
  
  # Print out zone-----
  out_qc_folder <- create_folder(out_qc_folder)
  out_plot_large <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-basic-large-plots.pdf"))
  out_plot_summary <- file.path(normalizePath(out_qc_folder), paste0(output_prefix,"-qc-basic-summary-plots.pdf"))
  
  .plotThePlots <- function(theplot, verbose = TRUE) {
    tryCatch(
      {
        print(theplot)
      },
      error=function(cond) {
        if(verbose) message("       - Plot cannot be printed")
        
      },
      warning=function(cond) {
        if(verbose) message(paste("       - Print density plots cause warnings"))
      })
  }
  
  if(printPDF){
    if(sn > 800){
      pdf(out_plot_large, width = 40, height = 8)
    }else if(sn <= 800 & sn > 200){
      pdf(out_plot_large, width = 22, height = 8)
    }else{
      pdf(out_plot_large, width = 14, height = 8)
    }
  }
  .plotThePlots(puid1)
  .plotThePlots(piseio)
  gridExtra::grid.arrange(piseso, puid2, ncol = 1, heights = c(2, 1))
  .plotThePlots(p_na_peprii)
  if(printPDF) garbage <- dev.off()
  
  if(printPDF) pdf(out_plot_summary, width = 12, height = 6)
  .plotThePlots(ptns)
  .plotThePlots(den1)
  .plotThePlots(den2)
  .plotThePlots(den3)
  if(printPDF) garbage <- dev.off()
}


