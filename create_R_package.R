#!/usr/bin/env Rscript

## create HCGB.IGTP.DAnalysis
library(devtools)
devtools::install_github("klutometis/roxygen")
setwd("../")
getwd()
devtools::create("HCGB.IGTP.DAnalysis")
