% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{df.factorizer}
\alias{df.factorizer}
\title{Converts columns in factors}
\usage{
df.factorizer(given.df, col_names.given, mode = "factor")
}
\arguments{
\item{given.df}{Dataframe}

\item{col_names.given}{List of columns to convert}

\item{mode}{Type of conversion to do: factor, as.factor, as.numeric = as.numeric(as.character("") )}
}
\description{
This functions loads a given dataframe and returns for the given set of columns, columns converted as factors, numeric, characters as specified.
}
