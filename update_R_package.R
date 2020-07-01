#!/usr/bin/env Rscript

## update XICRA.stats package
library(devtools)
devtools::install_github("klutometis/roxygen")

setwd("./")
getwd()
devtools::document()

