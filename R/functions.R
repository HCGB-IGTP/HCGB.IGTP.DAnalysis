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

#' Convert hg19 <-> hg38 coordinates
#' 
#' Using rtracklayer and GenomicRanges, this function converts hg19 <-> hg38 coordinates
#' @param df2convert Coordinate Dataframe provided to convert to GenomicRanges
#' @param folderWithLiftoverInfo Folder containing liftover chain information retrieved from UCSC website
#' @param hg19ToHg38 Boolean to transfrom hg19 to hg38
#' @param hg38ToHg19 Boolean to transfrom hg38 to hg19
#' @export

litfover_func <- function(df2convert, folderWithLiftoverInfo, 
                          hg19ToHg38=FALSE, hg38ToHg19=FALSE) {
  library(rtracklayer)
  library(GenomicRanges)
  
  ## create GRanges
  df2convert.GR <- GenomicRanges::makeGRangesFromDataFrame(df = df2convert)
  
  ## Check info:
  if (hg19ToHg38) {
    chain=import.chain(file.path(folderWithLiftoverInfo, "hg19ToHg38.over.chain")) # file downloaded from UCSC
  } else if (hg38ToHg19) {
    chain=import.chain(file.path(folderWithLiftoverInfo, "hg38ToHg19.over.chain")) # file downloaded from UCSC
  } else {
    print("ERROR: No option provided")
  }
  
  ## convert
  newLocs=liftOver(x = df2convert.GR, chain = chain)
  newLocsDF=data.frame(newLocs)
  
  ## check if there are some duplicates
  #rownames(newLocsDF) <- newLocsDF$group_name
  n_occur <- data.frame(table(newLocsDF$group_name))
  #df2convert.GR[n_occur[n_occur$Freq > 1,]$Var1,]
  #subset(newLocsDF, group_name %in% n_occur[n_occur$Freq > 1,]$Var1)
  
  rownames(newLocsDF) <- paste(newLocsDF$group_name, ave(newLocsDF$group_name, newLocsDF$group_name, 
                                                         FUN=function(i) seq(length(i))), sep='.')
  #subset(newLocsDF, group_name %in% n_occur[n_occur$Freq > 1,]$Var1)
  
  ##
  if (hg19ToHg38) {
    df2convert = data.frame(df2convert,"start.hg19"=df2convert$GeneStart)
    df2convert = data.frame(df2convert,"end.hg19"=df2convert$GeneEnd)
    df2convert[newLocsDF$group_name,"start.hg38"]=newLocsDF[,"start"]
    df2convert[newLocsDF$group_name,"end.hg38"]=newLocsDF[,"end"]
  } else if (hg38ToHg19) {
    df2convert = data.frame(df2convert,"start.hg38"=df2convert$GeneStart)
    df2convert = data.frame(df2convert,"end.hg38"=df2convert$GeneEnd)
    df2convert[newLocsDF$group_name,"start.hg19"]=newLocsDF[,"start"]
    df2convert[newLocsDF$group_name,"end.hg19"]=newLocsDF[,"end"]
  }
  
  
  return(df2convert)
}

#' Loads R data into variable
#' 
#' This functions loads a given RData object in a temporal environment and returns it
#' @param file2load Absolute path to save RData file
#' @export
loader <- function(file2load) {
  ## new_name <- loader(file2load = path_to_RData ))
  load(file = file2load,
       name2load <- new.env())
  return(name2load)
}

#' Saves R sessionInfo data into file
#' 
#' This functions saves the configuration of the session into a txt file
#' @param dir.given Absolute path to save file
#' @export
sessionInfo_write <- function(dir.given) {
  writeLines(capture.output(sessionInfo()), file.path(dir.given, "sessionInfo.txt"))  
}

#' Convert Spanish date: 14/5/2022 into English american format: 2022-5-14
#' 
#' @param date_str String date to convert
#' @export
convertDate <- function(date_str) {
  ## convert as data from Spanish Date to English Format: YYYYmmdd
  ## 15-04-2022 -> 2022-04-15
  return(str_split(date_str, "/") %>% unlist() %>% rev() %>% paste(collapse = '-'))
  
  ## Example: sapply(example_df$Date, convertDate)
}

#' Converts columns in factors
#' 
#' This functions loads a given dataframe and returns for the given set of columns, columns converted as factors, numeric, characters as specified.
#' @param given.df Dataframe
#' @param col_names.given List of columns to convert
#' @param mode Type of conversion to do: factor, as.factor, as.numeric = as.numeric(as.character("") )
#' @export
df.factorizer <- function(given.df, col_names.given, mode="factor"){
  
  if (mode=="factor") {
    # to do it for some names in a vector named 'col_names'
    given.df[col_names.given] <- lapply(given.df[col_names.given] , factor)
    
  } else if (mode=="as.factor") {
    # to do it for some names in a vector named 'col_names'
    given.df[col_names.given] <- lapply(given.df[col_names.given] , as.factor)
  }  else if (mode=="as.numeric") {
    # to do it for some names in a vector named 'col_names'
    given.df[col_names.given] <- lapply(given.df[col_names.given] , as.character)
    given.df[col_names.given] <- lapply(given.df[col_names.given] , as.numeric)
  }
  return(given.df)
}

