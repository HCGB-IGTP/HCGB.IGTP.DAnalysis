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


#' Save plots in PNG format
#'
#' @param folder_path 
#' @param name_file 
#' @param plot_given 
#'
#' @export
#'
save_png <- function(folder_path, name_file, plot_given) {
  png(file=file.path(folder_path, paste0(name_file, ".png")),
      width=1200, height=700)
  print(plot_given)  
  dev.off()
  
}

#' Save plots in JPEG format
#'
#' @param folder_path 
#' @param name_file 
#' @param plot_given 
#'
#' @export
#'
save_jpeg <- function(folder_path, name_file, plot_given) {
  jpeg(filename = file.path(folder_path, paste0(name_file, ".jpeg")),
       quality = 100,  width=1200, height=700)
  print(plot_given)  
  dev.off()
  
}

#' Save plots in multiple formats
#'
#' @param folder_path 
#' @param name_file 
#' @param plot_given 
#'
#' @export
#'
save_plots_multiformat <- function (folder_path, name_file, plot_given) 
{
  print("Save plots in several formats")
  save_pdf(folder_path, name_file, plot_given = plot_given)
  save_jpeg(folder_path, name_file, plot_given = plot_given)
  save_png(folder_path, name_file, plot_given = plot_given)
}

#' Save multiple plots in PDF
#'
#' @param folder_path Folder to store results
#' @param name_file Name of the file
#' @param list_plots List of plots to save in pdf
#'
#' @export
#'
save_multi_pdf <- function(folder_path, name_file, list_plots) {
  pdf(file.path(folder_path, paste0(name_file, ".pdf")), paper = "A4r",  width = 35, height = 12)
  for (i in list_plots) {
    print(i)
  }
  dev.off()
}
