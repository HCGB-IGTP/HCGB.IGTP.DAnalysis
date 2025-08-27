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


create_col_palette2 <- function(columnGiven, palette_given = "Paired") {
  
  levels_given <- levels(as.factor(columnGiven))
  library(RColorBrewer)
  colfactors <- factor(as.factor(columnGiven), levels = levels_given)
  n_colors <- length(levels(colfactors))
  if (n_colors < 3) {
    n_colors = 3
  }
  list_colors <- brewer.pal(n = n_colors, palette_given)
  
  toReturn <- list(
    "color.vector" = list_colors[colfactors], 
    "named.vector" = setNames(levels_given, list_colors[1:length(levels_given)] )
  )
  
  return(toReturn)
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

#' Get the length of unique elements in a vector
#'
#' @param vect_ Vector of elements
#'
#' @export
#'
get_length_vect <- function(vect_) { vect_ %>% as.factor() %>% levels() %>% length() }

#' Get the unique elements in a vector
#'
#' @param vect_ Vector of elements
#'
#' @export
#'
get_vect <- function(vect_) { vect_ %>% as.factor() %>% levels() }

#' Get the frequency of elements in a vector
#'
#' @param vect_ Vector of elements
#'
#' @export
get_freq <- function(vect_) { 
  library(tidyverse)
  df_ <- vect_ %>% as.factor() %>% table() %>% as.data.frame() 
  colnames(df_)[1] <- "Category"
  df_['%'] <- (df_$Freq/sum(df_$Freq))*100
  print(df_)
}

#' Get the not included in
#'
#' @param x Vector of elements 1
#' @param y Vector of elements 2 
#'
#' @export
'%!in%' <- function(x,y)!( '%in%'(x,y) )

#' Get the last element of a vector
#'
#' @param x Vector
#'
#' @export
last <- function(x) { return( x[length(x)] ) }



