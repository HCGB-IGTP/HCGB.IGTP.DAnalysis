% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RnBeads_functions.R
\name{load_beta_values}
\alias{load_beta_values}
\title{Load Beta Values}
\usage{
load_beta_values(betas_table_file, allow_missing_p = 0.95)
}
\arguments{
\item{betas_table_file}{Absolute path to CSV file}

\item{allow_missing_p}{Default=0.95. Allows only 5\% missing data}
}
\description{
This functions Reads into a dataframe the Beta values CSV file from RnBeads Track and table export module.
Removes missing rows using delete_na()
}
