
#' Generate Basic LAB QC Plots
#'
#' This function creates a set of four QC plots for an LAB assay. The plots include:
#' \itemize{
#'   \item An overall value distribution histogram with a density overlay and missing value percentage.
#'   \item A boxplot comparing value distributions across sample types.
#'   \item A scatter plot with a loess smooth showing the trend of values over injection or sample order.
#'   \item A bar plot summarizing the percentage of missing data by sample type.
#' }
#'
#' @param results A data frame containing assay results. Expected to have one row per analyte, with
#'   columns such as \code{analyte_name} and other metadata.
#' @param results_long A long-format data frame containing sample-level data. Expected columns include:
#'   \code{sample_id}, \code{sample_type}, \code{value}, \code{sample_order}, and \code{sample_id_ordered}.
#' @param out_qc_folder Character. The folder where QC plot PDFs will be saved. If \code{NULL},
#'   the function will create a default folder.
#' @param output_prefix Character. A prefix to be used in the names of output PDF files.
#' @param printPDF Logical. If \code{TRUE}, the function outputs the plots to PDF files.
#' @param verbose Logical. If \code{TRUE}, the function prints progress messages to the console.
#'
#' @return Invisibly returns a \code{gridExtra::grid.arrange} object containing the arranged plots.
#'   The primary output is the saved PDF file(s) if \code{printPDF = TRUE}.
#'
#' @details The function uses \code{ggplot2} syntax to generate clean and informative plots,
#'   avoiding clutter even with large numbers of samples (e.g., > 1400 samples). It calculates an
#'   overall missing value percentage for the measured values, and provides visual summaries of the
#'   data distribution and missing data prevalence.
#' @export
plot_basic_lab_qc <- function(results, 
                              results_long,
                              out_qc_folder = NULL,
                              output_prefix,
                              printPDF = TRUE,
                              verbose = TRUE) {
  
  sample_id = sample_type = value = sample_order = sample_id_ordered = NULL
  
  # Calculate overall missing percentage for 'value'
  total_samples <- base::nrow(results_long)
  missing_count <- sum(is.na(results_long$value))
  overall_missing <- round((missing_count / total_samples) * 100, 2)
  
  # Plot 1: Overall Value Distribution (Histogram with Density Overlay)
  p1 <- ggplot(results_long, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, 
                            fill = "skyblue", color = "black", na.rm = TRUE) +
    geom_density(color = "darkblue", linewidth = 1, na.rm = TRUE) +
    theme_minimal() +
    labs(title = "Overall Value Distribution",
                  subtitle = paste("Overall missing:", overall_missing, "%"),
                  x = "Value", y = "Density")
  
  # Plot 2: Value Distribution by Sample Type (Boxplot)
  p2 <- ggplot(results_long, aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(na.rm = TRUE) +
    theme_classic() +
    labs(title = "Value Distribution by Sample Type",
                  x = "Sample Type", y = "Value") +
    theme(legend.position = "none")
  
  # Plot 3: Value Trend Over Injection Order
  if (all(!is.na(as.numeric(as.character(results_long$sample_order))))) {
    results_long$sample_order_num <- as.numeric(as.character(results_long$sample_order))
    p3 <- ggplot(results_long, aes(x = sample_order_num, y = value, color = sample_type)) +
      geom_point(alpha = 0.5, na.rm = TRUE) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, na.rm = TRUE)
      theme_minimal() +
      labs(title = "Value Trend Over Injection Order",
                    x = "Injection Order", y = "Value")
  } else {
    p3 <- ggplot(results_long, aes(x = sample_id_ordered, y = value, color = sample_type)) +
      geom_point(alpha = 0.5, na.rm = TRUE) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, na.rm = TRUE)
      theme_minimal() +
      labs(title = "Value Trend Over Sample Order",
                    x = "Sample Order", y = "Value") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Plot 4: Missing Data Summary by Sample Type
  missing_summary <- results_long %>%
    dplyr::group_by(sample_type) %>%
    dplyr::summarise(total = dplyr::n(),
                     missing = sum(is.na(value))) %>%
    dplyr::mutate(missing_pct = round((missing / total) * 100, 2))
  
  p4 <- ggplot(missing_summary, aes(x = sample_type, y = missing_pct, fill = sample_type)) +
    geom_bar(stat = "identity", na.rm = TRUE) +
    theme_minimal() +
    labs(title = "Missing Data Percentage by Sample Type",
                  x = "Sample Type", y = "Missing (%)") +
    theme(legend.position = "none")

  if(verbose) message("  + (p) Plot QC plots")
  grid_plots <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  # Save to PDF if requested
  if (printPDF) {
    out_qc_folder <- create_folder(out_qc_folder)
    out_file <- file.path(normalizePath(out_qc_folder), paste0(output_prefix, "-lab-qc-plots.pdf"))
    
    pdf(out_file, width = 14, height = 10)
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
    dev.off()
  } else {
    print(grid_plots)
  }
}
