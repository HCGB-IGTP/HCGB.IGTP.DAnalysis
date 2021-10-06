#' Plots UpSetR
#' 
#' Creates plot with multiple interesection between sets using UpSetR package
#' @param data_set A dataframe containing as columns de sets, as rows the items and values are 0/1 for absence/presence. You can also provide create_upset_data(list_files)
#' @param sets Names of the set of interest to include in the comparison
#' @param y.label Name to include in the Y axis. Default: Items shared
#' @param x.label Name to include in the X axis. Default: Items/group
#' @return Returns plot and dataset used
#' @export
create_upset_plot <- function(data_set, sets, y.label="Items shared", x.label="Items/group") {
  library(UpSetR)
  p <- upset(data_set, sets = sets,
             mainbar.y.label = y.label,    # items shared
             sets.x.label = x.label,       # items/group
             order.by = "freq", sets.bar.color = "darkblue",
             point.size = 4,    matrix.color = "Red",
             nintersects = 100, text.scale = 1.5,  keep.order = TRUE) ## order set group provided
  print(p)
  
  list2return <- list(
    "plot" = p,
    "dataset" = data_set
  )
  return(list2return)
}

#' Create UpSetR plot dataframe
#' 
#' Create dataframe with presence/abscence from list to plot interesection between sets using UpSetR package
#' @param list_files A list of names list of characters. 
#' @param sets Names of the set of interest to include in the comparison
#' @param y.label Name to include in the Y axis. Default: Items shared
#' @param x.label Name to include in the X axis. Default: Items/group
#' @export
create_upset_data <- function(list_files) {
  list_files.ids <- unique(unlist(list_files))
  df_dat <- as.data.frame(sapply(list_files, function(x) table(factor(x, levels=list_files.ids))))
  return(df_dat)
}


#' Get files for UpSetR plot
#'  
#' Get information of interest to produce UpSetR plots from a folder containing files
#' @param data_dir Absolute path to folder containing files with e.g. comp1-comp2_down.whatever.txt 
#' @param pattern2search Pattern to include in regex to search file names. e.g. down, complex.isomir, etc
#' @return List of lists containing data for each set
#' @export
get_data_upset <- function(data_dir, pattern2search) {
  list_files <- list.files(data_dir, pattern2search)
  listdata = list()
  for (item in list_files) {
    print(item)
    list_produced <- unlist(strsplit(item, split="\\.")) ## comp1-comp2_down.whatever.txt 
    str_name = list_produced[1] #comp1-comp2_down
    listdata[[str_name]] = readLines(file.path(data_dir, item))
  }
  ## 
  return(listdata)
}

#' Creates UpSetR plot
#' 
#' Get information of interest to produce UpSetR plots from a folder containing files and return plot generated
#' @param data_dir Absolute path to folder containing files with e.g. comp1-comp2_down.whatever.txt 
#' @param pattern2search Pattern to include in regex to search file names. e.g. down, complex.isomir, etc
#' @return List of plot and data containing data for each set
#' @export
create_upset <- function(data_dir, pattern2search) {
  
  ## get data   
  list_files <- get_data_upset(data_dir, pattern2search)
  
  ## create plot
  returnData <- create_upset_plot(data_set = create_upset_data(list_files), 
                                  sets = names(list_files))
  
  data2return <- list(
    "upset_plot" = returnData,
    "listFiles" = list_files
  )
  
  return(data2return)
}