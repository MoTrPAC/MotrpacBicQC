% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolomics_qc.R
\name{load_metabolomics_batch}
\alias{load_metabolomics_batch}
\title{Load metabolomics batch}
\usage{
load_metabolomics_batch(input_results_folder, cas, verbose = TRUE)
}
\arguments{
\item{input_results_folder}{(char) Path to the PROCESSED_YYYYMMDD folder}

\item{cas}{(char) Chemical Analytical Site code (e.g "umichigan")}

\item{verbose}{(logical) \code{TRUE} (default) shows messages}
}
\value{
(list of data.frames) List of all the data frames
}
\description{
Open, check, and return all metabolomics files
}
\examples{
\dontrun{
here <- load_metabolomics_batch(input_results_folder = "/path/to/PROCESSED_YYYYMMDD/", 
                                cas = "cassite")
}
}
