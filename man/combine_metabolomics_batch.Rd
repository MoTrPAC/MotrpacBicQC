% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolomics_qc.R
\name{combine_metabolomics_batch}
\alias{combine_metabolomics_batch}
\title{Combines all files from a Metabolomics batch}
\usage{
combine_metabolomics_batch(input_results_folder, cas, verbose = TRUE)
}
\arguments{
\item{input_results_folder}{(char) Path to the PROCESSED_YYYYMMDD folder}

\item{cas}{(char) Chemical Analytical Site code (e.g "umichigan")}

\item{verbose}{(logical) \code{TRUE} (default) shows messages}
}
\value{
(int) number of issues identified
}
\description{
Combines all the files from an untargeted assay submission
including metadata_sample, metadata_metabolites, results for both "NAMED" and "UNNAMED"
folders
}
\examples{
\dontrun{
all_datasets <- combine_metabolomics_batch(
                        input_results_folder = "/full/path/to/PROCESSED_YYYYMMDD/",
                        cas = "umichigan")
}
}
