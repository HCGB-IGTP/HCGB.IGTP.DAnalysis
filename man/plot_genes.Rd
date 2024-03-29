% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RnBeads_functions.R
\name{plot_genes}
\alias{plot_genes}
\title{RnBeads Locus plotter dataframe}
\usage{
plot_genes(
  data_frame_given,
  name_given,
  RnBead_object,
  sample.grouping,
  out_folder
)
}
\arguments{
\item{data_frame_given}{Dataframe with information from differential_methylation_data/diffMethTable_region_genes...}

\item{name_given}{name to use}

\item{sample.grouping}{rnb.sample.groups output grouping samples}

\item{out_folder}{Folder to store pdfs}

\item{RnBeads_object}{RnBeads object}
}
\description{
This functions calls locus_plotter (rnb.plot.locus.profile caller) for each entry in a dataframe
}
