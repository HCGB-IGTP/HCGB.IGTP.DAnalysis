#' Get numbers rows fulfilling cutoff of filt_length
#'
#' @param i 
#' @param filt_length 
#'
#' @export
get_dim.filt <- function(i, filt_length=1e4) { 
  
  i.df <- as.data.frame(i)
  i.subset <- subset(i.df, length > filt_length)
  return(dim(i.subset)[1] )
}

#' Get categories of gene length regions
#'
#' @param df_given Dataframe to use containing a width or length column
#'
#' @export
get_length_cat <- function(df_given, col_name="width") {
  
  length2use = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
  names2use = c("1kbp", "10kbp", "100kbp", "1Mbp", "10Mbp", "100Mbp")
  
  df_given['length_cat'] <- lapply( df_given[[ col_name ]], function(x) {
        paste(names2use[x > length2use], collapse = "-")
  }) %>% unlist()
  
  
  ## rename categories
  ## [1] ""                             "1kbp"        
  ## [3] "1kbp-10kbp"                   "1kbp-10kbp-100kbp"           
  ## [5] "1kbp-10kbp-100kbp-1Mbp"       "1kbp-10kbp-100kbp-1Mbp-10Mbp"
  ## 
  df_given <- df_given %>% mutate(length_cat=case_when( length_cat == "" ~ "<10kbp",
                                                        length_cat == "1kbp" ~ "1-10kbp",
                                                        length_cat == "1kbp-10kbp" ~ "10-100kbp",
                                                        length_cat == "1kbp-10kbp-100kbp" ~ "0.1-1Mbp",
                                                        length_cat == "1kbp-10kbp-100kbp-1Mbp" ~ "1-10Mbp",
                                                        length_cat == "1kbp-10kbp-100kbp-1Mbp-10Mbp" ~ "10-100Mbp"
  ))
  
  ## Order categories
  df_given$length_cat <- factor(df_given$length_cat, levels = c("<10kbp", "1-10kbp", "10-100kbp",
                                                                "0.1-1Mbp", "1-10Mbp", "10-100Mbp"))
  
  ## print frequencies
  df_given$length_cat %>% get_freq()
  
  
  df_given
  
}

#' Get sum of bp between min and max length
#'
#' @param df_given Dataframe to use containing a length column
#'
#' @export
get_sum_bp <- function (i, min_length = 100, max_length=1000) {
  i.df <- as.data.frame(i)
  i.subset <- subset(i.df, length > min_length & length <= max_length)
  return(sum(i.subset$length))
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

#' Get rows of dataframe
#'
#' @param obj_given Dataframe, matrix or 
#'
#' @export
#'
get_rows <- function(obj_given) { dim(as.data.frame(obj_given, row.names = NULL))[1] }


#' Get cols of dataframe
#'
#' @param obj_given Dataframe, matrix or 
#'
#' @export
get_cols <- function(obj_given) { dim(as.data.frame(obj_given, row.names = NULL))[2] }

#' Remove rownames
#'
#' @param df Dataframe
#'
#' @export
remove_rownames <- function(df) {
  df <- as.data.frame(df)
  df['sample'] <- rownames(df)
  rownames(df) <- NULL
  df
}

#' Get the last element of a vector
#'
#' @param df Dataframe
#' @param col_name Column name to use as new rownames
#'
#' @export
add_rownames <- function(df, col_name) {
  df <- as.data.frame(df)
  rownames(df) <- df[[ col_name ]]
  df
}