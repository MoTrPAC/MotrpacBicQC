
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot Basic Metabolomics QC charts
#'
#' @description Plot intensity distributions, number of unique ids, NA values 
#' and RT/MZ density maps (only untargeted)
#' @param results (df) metabolomics results (both named and unnamed already merged, if untargeted)
#' @param results_long (df) metabolomics results (both named and unnamed already merged, if untargeted, in long format)
#' @param metametab (df) metadata samples metadata (only untargeted sites)
#' @param out_qc_folder (char) output qc folder (it creates the folder if it doesn't exist)
#' @param output_prefix (char) prefix for the file name output (pdf file)
#' @param untargeted (logical) `TRUE` if the dataset is untargeted (named + unnamed metabolites)
#' @param printPDF (logical) `TRUE` (default) prints pdf file
#' @param verbose (logical) `TRUE` (default) shows messages
#' @export
plot_basic_metabolomics_qc <- function(results, 
                                       results_long,
                                       metametab = NULL,
                                       out_qc_folder = NULL,
                                       output_prefix,
                                       printPDF = TRUE,
                                       untargeted = TRUE,
                                       verbose = TRUE){

  metabolite_name = id_type = sample_id = sample_order = intensity = sample_type = sum_quant = mz = rt = ..density.. = NULL
  
  if(verbose) message("   + (+) QC PLOTS ------------------")
  
  if(verbose) message("       - (p) Plot intensity distributions")
  
  # Get the total number of samples to customize the plots. If larger than 200, 
  # prepare for large plots
  
  sn <- length(unique(results_long$sample_id))
  # Set a limit for which remove labels 
  sn_limit <- 200
  
  if(!is.null(metametab)){
    pmzrt <- ggplot(metametab, aes(x=mz, y=rt) ) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      viridis::scale_fill_viridis() +
      facet_wrap(~id_type) +
      theme(
        legend.position='none'
      ) +
      labs(title = "MZ/RT density map", 
           subtitle = paste(output_prefix))
  }
  
  # Density plots------ 
  
  if(verbose) message("       - (p) Density distributions")
  mu <- results_long %>%
    group_by(sample_id) %>%
    reframe(grp.mean = mean(intensity))
  
  den1 <- ggplot(data = results_long, aes(x = log2(intensity), color = sample_type)) + 
    geom_density(na.rm = TRUE) + 
    theme_light() +
    labs(title = "Density distribution by sample type", 
         subtitle = paste(output_prefix)) 
  
  den2 <- ggplot(data = results_long, aes(x = log2(intensity), 
                                          color = sample_id)) +
    geom_density(na.rm = TRUE) + 
    theme_light() +
    theme(legend.position="none") +
    labs(title = "Density distribution by sample id", 
         subtitle = paste(output_prefix))
  
  den3 <- ggplot(data = results_long, aes(x = log2(intensity), 
                                          color = sample_id)) +
    geom_density(na.rm = TRUE) + 
    theme_minimal() +
    theme(legend.position="none") + 
    facet_wrap(~sample_type) +
    labs(title = "Density distributions", 
         subtitle = paste(output_prefix)) 
  
  # Plot: sum of intensities------
  if(verbose) message("       - (p) Plot sum of intensity/concentration values")
  sum_int <- results_long %>%
    group_by(sample_id, sample_type, sample_order) %>%
    reframe(sum_quant = sum(intensity))
  
  psumint <- ggplot(sum_int, aes(x = reorder(sample_id, sample_order), y = sum_quant, fill = sample_type)) +
    geom_bar(stat = "identity") + theme_classic() +
    theme(axis.text.y = element_text(angle = 90), legend.position="top")
  if(untargeted){
    psumint <- psumint + labs(
      title = "Sum of Intensities by sample", 
      subtitle = paste(output_prefix), 
      x = "Samples",
      y = "Sum of intensity"
    )
  }else{
    psumint <- psumint + labs(
      title = "Sum of concentrations by sample", 
      subtitle = paste(output_prefix), 
      x = "Samples",
      y = "Sum of Concentrations"
    )
  }
  if(sn > sn_limit){
    psumint <- psumint + theme(
      axis.text.x = element_blank())
  }else{
    psumint <- psumint + theme(
      axis.text.x = element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5,
                                 size = 6))
  }
  
  # Boxplot distribution, as part of a grid----
  piseso <- ggplot2::ggplot(results_long,
                            aes(
                              x = reorder(sample_id, sample_order),
                              y = log2(intensity),
                              fill = sample_type)) +
    geom_boxplot(na.rm = TRUE) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 90),
      legend.position="top") +
    if(untargeted){
      labs(
        title = "Intensity distribution / Unique IDs", 
        subtitle = paste(output_prefix), 
        y = "log2(intensity)"
        )
    }else{
      labs(
        title = "Concentration distribution / Unique IDs", 
        subtitle = paste(output_prefix), 
        y = "log2(Concentration)"
      )
    }
  
  #Boxplots: standalone----
  piseio <- 
    ggplot2::ggplot(results_long,
                    aes(x = reorder(sample_id, sample_order),
                        y = log2(intensity),
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
    if(untargeted){
      piseio <- piseio + labs(title = "Sample Intensity distribution",
                              y = "log2(intensity)") 
    }else{
      piseio <- piseio + labs(title = "Sample Concentration distribution",
                              y = "log2(Concentration)") 
    }
  
  # Plot number fo IDs by identification------
  if(verbose) message("       - (p) Plot ID counts")
  uid <- results_long %>%
    group_by(across(all_of(c("metabolite_name", "sample_id", "sample_type", "sample_order", "id_type")))) %>%
    reframe(total_intensity = intensity)
  uid2 <- uid[which(!is.na(uid$total_intensity)),]
  uid3 <- unique(uid2[c("metabolite_name", "sample_id", "sample_type", "sample_order", "id_type")]) %>%
    count(sample_id, sample_type, sample_order, id_type)
    
  
  puid1 <- ggplot(uid3, 
                  aes(x = reorder(sample_id, sample_order), 
                      y = n, 
                      fill = sample_type)) +
    geom_bar(stat = "identity",
             position="stack",
             na.rm = TRUE) +
    theme_classic() +
    labs(title = "Unique IDs in samples", 
         subtitle = paste(output_prefix),
         caption = "Scales are different",
         x = "Sample ID", y = "Number of identifications") +
    facet_grid(rows = vars(id_type)) +
    theme(legend.position = "top")
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
  uid4 <- unique(uid2[c("metabolite_name", "sample_id", "sample_type", "sample_order")]) %>%
    count(sample_id, sample_type, sample_order)
  
  puid2 <- ggplot(uid4, 
                  aes(x = reorder(sample_id, sample_order), 
                      y = n, 
                      fill = sample_type)) +
  geom_bar(stat = "identity",
           na.rm = TRUE) +
    theme_classic() +
    labs(caption = "Sort by injection order",
         x = "Sample ID", 
         y = "") +
    theme(legend.position = "none")
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
  tns <- unique(uid2[c("sample_id", "sample_type")]) %>% group_by(sample_type) %>% count(sample_type)
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
  
  
  ppids <- ggplot(results, aes(x = as.factor(id_type), fill = as.factor(id_type))) +
    geom_bar(aes(y = (..count..)/sum(..count..))) +
    geom_text(aes(y = ((..count..)/sum(..count..)), 
                  label = scales::percent((..count..)/sum(..count..))), 
              stat = "count", vjust = -0.25) +
    scale_y_continuous(labels = percent) +
    labs(title = "Proportion of Features Identified (named/unnamed)", 
         subtitle = paste(output_prefix), 
         y = "", x = "") +
    theme_light() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18)) + 
    scale_fill_brewer(palette = "Dark2")
  
  pnids <- ggplot(results, aes(x = id_type, fill = id_type)) +
    geom_bar() +
    geom_text(stat = 'count',aes(label =..count.., vjust = -0.2)) +
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
  .plotThePlots(psumint)
  .plotThePlots(puid1)
  .plotThePlots(piseio)
  gridExtra::grid.arrange(piseso, puid2, ncol = 1, heights = c(2, 1))
  .plotThePlots(p_na_peprii)
  if(printPDF) garbage <- dev.off()
  
  if(printPDF) pdf(out_plot_summary, width = 12, height = 6)
  .plotThePlots(ptns)
  gridExtra::grid.arrange(ppids, pnids, ncol = 2)
  .plotThePlots(den1)
  .plotThePlots(den2)
  .plotThePlots(den3)
  if(!is.null(metametab)){.plotThePlots(pmzrt)}
  if(printPDF) garbage <- dev.off()
}


