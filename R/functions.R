#' Assigns color codes for samples 
#' 
#' Create palette color codes for each sample given using RColorBrewer given a variable
#' @param columnGiven Column dataframe given to use as class for each color
#' @param levels_given Orderered or unordered levels to use
#' @param palette_given Name of the palette to use: Paired as default.
#' @export
create_col_palette <- function(columnGiven, levels_given, palette_given="Paired") {
  
  library(RColorBrewer)
  ## reorder levels:
  colfactors <- factor(as.factor(columnGiven), levels=levels_given)
  
  ## set number of levels to use
  n_colors <- length(levels(colfactors))
  if (n_colors<3) { n_colors=3 } ## Fix: Use at least 3
  
  ## create color code
  list_colors <- brewer.pal(n=n_colors, palette_given)
  
  return(list_colors[colfactors])
}

#' Save plot in pdf
#' 
#' Save a plot in a pdf in the folder provided.
#' @param folder_path Absolute path to save pdf
#' @param name_file PDF file name, do not include extension
#' @param plot_given Plot object
#' @export
save_pdf <- function(folder_path, name_file, plot_given) {
  pdf(file.path(folder_path, paste0(name_file, ".pdf")), paper = "A4r", width = 35, height = 12)
  print(plot_given)
  dev.off()
}

#' Remove NAs in dataframe
#' 
#' Filters out dataframe according the amount of NAs allowed
#' @param DF Dataframe provide
#' @param n Number of NAs allowed for each row.
#' @export
delete_na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}


