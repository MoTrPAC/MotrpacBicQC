% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots_olink.R
\name{plot_basic_olink_qc}
\alias{plot_basic_olink_qc}
\title{Plot Basic OLINK QC charts}
\usage{
plot_basic_olink_qc(
  results,
  results_long,
  out_qc_folder = NULL,
  output_prefix,
  printPDF = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{results}{(df) olink results merged with proteins metadata}

\item{results_long}{(df) olink results in long format and merged with
sample metadata}

\item{out_qc_folder}{(char) output qc folder (it creates the folder if it doesn't exist)}

\item{output_prefix}{(char) prefix for the file name output (pdf file)}

\item{printPDF}{(logical) \code{TRUE} (default) prints pdf file}

\item{verbose}{(logical) \code{TRUE} (default) shows messages}
}
\description{
Plot value distributions, number of unique ids, NA values
}
