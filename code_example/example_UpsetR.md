
############################
## UpSetR
############################

## install UpSetR if not installed and loaded it
install.packages("UpSetR")
library(UpSetR)
############################

############################
## HCGB.IGTP.DAnalysis
############################
## I have generated some function to speed the use of UpsetR
## they are hosted in github here:
## https://github.com/HCGB-IGTP/HCGB.IGTP.DAnalysis
## and specifically here:
## https://github.com/HCGB-IGTP/HCGB.IGTP.DAnalysis/blob/master/R/UpSetR_functions.R

## Try to install HCGB.IGTP.DAnalysis or just copy paste functions
## install devtools if not installed
install.packages("devtools")

## install HCGB R package
devtools::install_github("HCGB-IGTP/HCGB.IGTP.Danalysis")

## load it
library(HCGB.IGTP.Danalysis)


## If it fails, here are the functions necessary

##----------------------------
## functions
##----------------------------
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
  
  print(list_files)
  
  ## create plot
  returnData <- create_upset_plot(data_set = create_upset_data(list_files), 
                                  sets = names(list_files))
  
  data2return <- list(
    "upset_plot" = returnData,
    "listFiles" = list_files
  )
  
  return(data2return)
}

#' Get UpSetR unique intersect
#' 
#' Get information of the unique items in a given intersect
#' @param data_set Matrix containing information for each set and item. Created from: create_upset_data() and/or returned by create_upset()
#' @param set2test Set name to retrieve unique single items. It can be one or two items. e.g. male-female_down, c("example_down", "example-cond2_down")
#' @return List of items unique for the set of interest
#' @export
get_single_UpSet_plot <- function(data_set, set2test) {
  data_df <- data.frame()
  if (length(set2test)==2) {
    data_df <- data_set[ data_set[[set2test[1]]]==1 & data_set[[set2test[1]]] == data_set[[set2test[2]]] & rowSums(data_set)==2,]  
  } else if (length(set2test)==1) {
    data_df <- data_set[ data_set[[set2test]]==1 & rowSums(data_set)==1,]  
  } else {
    print("ERROR: Option not available")
  }
  return(rownames(data_df))
}
##----------------------------
############################

############################
## Example to create UpsetR plot
############################

##----------------------------
## UpsetR from list of items
##----------------------------

## Create example set of lists
list.1 <- c("a", "b", "c", "d", "e", "f", "g")
list.2 <- c("a", "b", "c", "d", "e", "g")
list.3 <- c("a", "b", "c", "h", "m")
list.4 <- c("a", "b", "c", "g", "i", "j", "k", "l")

list.of.list <- list(
  "list.1" = list.1,
  "list.2" = list.2,
  "list.3" = list.3,
  "list.4" = list.4
)


## create a matrix/dataframe of presence/absence (0/1) for each
## list and each item: create_upset_data
df.presence <- create_upset_data(list.of.list)
df.presence

## Create plot: create_upset_plot
## Use df.presence generated and names for each set
create_upset_plot(data_set = df.presence, sets = names(df.presence))
##----------------------------

##----------------------------
## UpsetR from folder containing files with IDs, genes, names, whatever
##----------------------------
## Create plot from folder: read files, inherit names from file.names and create upset plot
## use function: create_upset

##----------------------------
## Save example list generated before
##----------------------------
test_upset.folder <- "test_UpsetR"
dir.create(test_upset.folder)

write.table(list.1, file = file.path(test_upset.folder, "list-1.example.txt"), 
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(list.2, file = file.path(test_upset.folder, "list-2.example.txt"), 
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(list.3, file = file.path(test_upset.folder, "list-3.example.txt"), 
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(list.4, file = file.path(test_upset.folder, "list-4.example.txt"), 
            row.names = FALSE, quote = FALSE, col.names = FALSE)
##----------------------------

##----------------------------
## Create:
# get files with given pattern, create dataframe, create matrix and create plot
upset_generated <- create_upset(data_dir = test_upset.folder, 
             pattern2search = ".example.txt")

upset_generated$upset_plot$dataset
upset_generated$upset_plot$plot
upset_generated$listFiles$`list-1`
##----------------------------
  