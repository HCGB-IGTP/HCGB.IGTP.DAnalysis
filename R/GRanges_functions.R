
#' Get dimensions for a GRanges object
#'
#' @param granges_given 
#'
#' @export
#'
dim_GRanges <- function(granges_given) {
  dim(as.data.frame(granges_given))
}

#' Get column names for a GRanges object
#'
#' @param granges_given 
#'
#' @export
#'
colnames_GRanges <- function(granges_given) {
  colnames(as.data.frame(granges_given))
}

#' Convert hg19 <-> hg38 coordinates
#' 
#' Using rtracklayer and GenomicRanges, this function converts hg19 <-> hg38 coordinates
#' @param df2convert GenomicRanges or coordinate Dataframe provided to convert to GenomicRanges
#' @param folderWithLiftoverInfo Folder containing liftover chain information retrieved from UCSC website
#' @param hg19ToHg38 Boolean to transfrom hg19 to hg38
#' @param hg38ToHg19 Boolean to transfrom hg38 to hg19
#' @export
litfover_func <- function(folderWithLiftoverInfo, data2convert,
                          hg19ToHg38=FALSE, hg38ToHg19=FALSE) {
  library(rtracklayer)
  library(GenomicRanges)
  
  if (is.data.frame(data2convert)) {
    
    ## create GRanges
    df2convert.GR <- GenomicRanges::makeGRangesFromDataFrame(df = data2convert)
  } else {
    ## Assume is a GRanges object already
    df2convert.GR <- data2convert
  }
  
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
  
  newLocs.GR <- unlist(newLocs)
  return(newLocs.GR)
  
  # newLocsDF=data.frame(newLocs)
  # 
  # ## check if there are some duplicates
  # #rownames(newLocsDF) <- newLocsDF$group_name
  # n_occur <- data.frame(table(newLocsDF$group_name))
  # #df2convert.GR[n_occur[n_occur$Freq > 1,]$Var1,]
  # #subset(newLocsDF, group_name %in% n_occur[n_occur$Freq > 1,]$Var1)
  # 
  # rownames(newLocsDF) <- paste(newLocsDF$group_name, ave(newLocsDF$group_name, newLocsDF$group_name, 
  #                                                        FUN=function(i) seq(length(i))), sep='.')
  # #subset(newLocsDF, group_name %in% n_occur[n_occur$Freq > 1,]$Var1)
  # 
  # ##
  # if (hg19ToHg38) {
  #   df2convert = data.frame(df2convert,"start.hg19"=df2convert$GeneStart)
  #   df2convert = data.frame(df2convert,"end.hg19"=df2convert$GeneEnd)
  #   df2convert[newLocsDF$group_name,"start.hg38"]=newLocsDF[,"start"]
  #   df2convert[newLocsDF$group_name,"end.hg38"]=newLocsDF[,"end"]
  # } else if (hg38ToHg19) {
  #   df2convert = data.frame(df2convert,"start.hg38"=df2convert$GeneStart)
  #   df2convert = data.frame(df2convert,"end.hg38"=df2convert$GeneEnd)
  #   df2convert[newLocsDF$group_name,"start.hg19"]=newLocsDF[,"start"]
  #   df2convert[newLocsDF$group_name,"end.hg19"]=newLocsDF[,"end"]
  # }
  # 
  # 
  # return(df2convert)
}

#' Subset lisf ot GRanges objects
#'
#' @param list_of_GRanges 
#' @param sub_Granges GRanges to use to subset 
#' @param string_given String to use to subset 
#'
#' @export
#'
GRanges_subsetter <- function(list_of_GRanges, sub_Granges=NULL, string_given=NULL) {
  
  ## store information as a list for each element of list_of_GRanges
  data2return <- list()
  
  ###################################
  ## check & control if 
  ###################################
  if (is.null(names(list_of_GRanges))) {
    # No list provided, a single GRanges
    print("ERROR: No list of GRanges provided")
    return()
  }
  
  if (is.null(sub_Granges)) {
    if (is.null(string_given)) {
      print("ERROR: No option provided")
      return()
    }
    print("Option 2: Use information to filter as a string")
    print("--------------------------")
    print(string_given)
  } else {
    print("Option 1: Use information saved as a GRanges object")
  }
  ###################################
  
  library(plyranges)
  library(IRanges)
  
  ### use a criteria: seqnames, start & stop or whatever provided as a string
  for (gr_sub in names(list_of_GRanges)) {
    
    print("Subsetting: ")
    print(gr_sub)
    
    if (is.null(sub_Granges)) {
      
      ###------------------------
      ### Use plyranges
      ###------------------------
      ### Use a given string: e.g. "seqnames == '1', start >= 1, end <= 1e6"
      ### we will split by commas and filter as many times as necessary
      
      list_String2use <- stringr::str_split(string_given, ",") %>% unlist()
      gr.tmp <- list_of_GRanges[[gr_sub]]
      for (ex in list_String2use) {
        gr.tmp <- gr.tmp %>% plyranges::filter(!! rlang::parse_expr(ex))
      }
      data2return[[ gr_sub ]] <- gr.tmp
      
    } else {
      
      ###------------------------
      ## Use subsetByOverlaps from IRanges
      ###------------------------
      ## use a GRanges object
      data2return[[ gr_sub ]] <- IRanges::subsetByOverlaps(x = list_of_GRanges[[gr_sub]], ranges = sub_Granges)
    }
    
  }  
  
  ###################################
  ## create summary
  ###################################
  df_dimensions.here <- lapply(data2return, dim_GRanges) %>% as.data.frame() %>% t()
  colnames(df_dimensions.here) <- c('filtered_rows', 'cols')
  #print(df_dimensions.here)
  
  df_dimensions_orig.here <- lapply(list_of_GRanges, dim_GRanges) %>% as.data.frame() %>% t()
  colnames(df_dimensions_orig.here) <- c('original_rows', 'cols')
  #print(df_dimensions_orig.here)
  df_dimensions_merge <- merge(df_dimensions_orig.here, df_dimensions.here, by="row.names")
  df_dimensions_merge$cols.x <- NULL
  rownames(df_dimensions_merge) <- df_dimensions_merge$Row.names
  df_dimensions_merge$Row.names <- NULL
  print(df_dimensions_merge)
  
  data2return[['dimensions']] = df_dimensions_merge
  
  ## save input information 
  data2return[['input']][['string']] = string_given
  data2return[['input']][['GRanges']] = sub_Granges
  
  
  return(data2return)
}

#' Create matrix of region overlaps
#' 
#' This functions takes a genomicRanges object and using regioneR::numOverlaps creates a matrix of 
#' overlapping regions
#'
#' @param GR.list_given List of Genomic ranges objects per sample
#'
#' @return matrix
#' @export
create_matrix_overlaps <- function(GR.list_given, n_cores=4) {
  
  matrixinp = matrix(data=0, nrow=length(GR.list_given), ncol=length(GR.list_given)) 
  
  for(i in 1:length(GR.list_given)){ 
    matrixinp[i,] = unlist(
      create_overlaps(GR.list_given[[i]], 
                      GR.list_given, 
                      mc.cores = n_cores))
  }
  
  # print(matrixinp) 
  colnames(matrixinp) <- names(GR.list_given)
  rownames(matrixinp) <- names(GR.list_given)
  
  return(matrixinp)
}

#' Create overlaps for a GR vs. list of GRs
#'
#' @param GR_given 
#' @param GR.list_given 
#' @param mc.cores 
#'
#' @export
create_overlaps <- function(GR_given, GR.list_given, mc.cores=2) { 
  mclapply(GR.list_given, FUN=regioneR::numOverlaps, B=GR_given, mc.cores = mc.cores)  
}





