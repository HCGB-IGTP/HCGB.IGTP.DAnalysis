#!/usr/bin/env Rscript

## update XICRA.stats package
#install.packages('desc')
library(devtools)
devtools::install_github("klutometis/roxygen")

setwd("./")
getwd()
devtools::document()

