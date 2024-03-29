% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validations.R
\name{validate_na_empty}
\alias{validate_na_empty}
\title{Validate Column for NA and Empty Values}
\usage{
validate_na_empty(df, col_name, verbose = TRUE)
}
\arguments{
\item{df}{A data frame.}

\item{col_name}{A character string specifying the name of the column to check.}

\item{verbose}{A logical indicating whether to print informative messages. Default is TRUE.}
}
\value{
Number of issues
}
\description{
This function checks if a specified column in a data frame contains either NA or empty values.
}
\examples{
df <- data.frame(A = c("a", "", NA, "d"), B = 1:4)
validate_na_empty(df, "A")

}
