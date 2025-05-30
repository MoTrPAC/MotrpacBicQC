% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_lab.R
\name{plot_basic_lab_qc}
\alias{plot_basic_lab_qc}
\title{Generate Basic LAB QC Plots}
\usage{
plot_basic_lab_qc(
  results,
  results_long,
  out_qc_folder = NULL,
  output_prefix,
  printPDF = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{results}{A data frame containing assay results. Expected to have one row per analyte, with
columns such as \code{analyte_name} and other metadata.}

\item{results_long}{A long-format data frame containing sample-level data. Expected columns include:
\code{sample_id}, \code{sample_type}, \code{value}, \code{sample_order}, and \code{sample_id_ordered}.}

\item{out_qc_folder}{Character. The folder where QC plot PDFs will be saved. If \code{NULL},
the function will create a default folder.}

\item{output_prefix}{Character. A prefix to be used in the names of output PDF files.}

\item{printPDF}{Logical. If \code{TRUE}, the function outputs the plots to PDF files.}

\item{verbose}{Logical. If \code{TRUE}, the function prints progress messages to the console.}
}
\value{
Invisibly returns a \code{gridExtra::grid.arrange} object containing the arranged plots.
The primary output is the saved PDF file(s) if \code{printPDF = TRUE}.
}
\description{
This function creates a set of four QC plots for an LAB assay. The plots include:
\itemize{
\item An overall value distribution histogram with a density overlay and missing value percentage.
\item A boxplot comparing value distributions across sample types.
\item A scatter plot with a loess smooth showing the trend of values over injection or sample order.
\item A bar plot summarizing the percentage of missing data by sample type.
}
}
\details{
The function uses \code{ggplot2} syntax to generate clean and informative plots,
avoiding clutter even with large numbers of samples (e.g., > 1400 samples). It calculates an
overall missing value percentage for the measured values, and provides visual summaries of the
data distribution and missing data prevalence.
}
