

#' Create matrix of region overlaps (by sequence)
#' 
#' This functions takes a genomicRanges object and using regioneR::numOverlaps creates a matrix of 
#' overlapping regions by sequence.
#'
#' @param GR.list_given List of Genomic ranges objects per sample
#' @param list_of_chrs List of seq IDs to use, by default chromosomes 1-22.
#'
#' @export
#'
create_matrix_overlaps_by_seqnames <- function(GR.list_given, list_of_chrs=c(1:22)) {
  
  list_of_chrs_matrix <- list()
  for (chr in list_of_chrs) {
    matrixinp = matrix(data=0, nrow=length(GR.list_given), ncol=length(GR.list_given)) 
    
    # fill the elements with j values 
    # in a matrix 
    for(j in 1:length(GR.list_given)){ 
      for(i in 1:length(GR.list_given)){ 
        
        i_GR <- GR.list_given[[i]]
        j_GR <- GR.list_given[[j]]
        matrixinp[i,j] = regioneR::numOverlaps(i_GR[seqnames(i_GR)==chr], j_GR[seqnames(j_GR)==chr])  
      } 
    } 
    
    # print(matrixinp) 
    colnames(matrixinp) <- names(GR.list_given)
    rownames(matrixinp) <- names(GR.list_given)
    
    matrixinp_df <- as.data.frame(matrixinp)
    matrixinp_df['chr'] <- chr
    
    list_of_chrs_matrix[[ paste0('chr_', chr) ]] <- matrixinp_df
  } 
  
  return(list_of_chrs_matrix)
}

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
#' @param verbose Print comments or not
#'
#' @export
#'
GRanges_subsetter <- function(list_of_GRanges, sub_Granges=NULL, string_given=NULL, verbose=TRUE) {
  
  ## store information as a list for each element of list_of_GRanges
  data2return <- list()
  
  ###################################
  ## check & control if 
  ###################################
  if (verbose) {
    
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
  }
    ###################################
  
  library(plyranges)
  library(IRanges)
  
  ### use a criteria: seqnames, start & stop or whatever provided as a string
  for (gr_sub in names(list_of_GRanges)) {
    
    if (verbose) {
      print("Subsetting: ")
      print(gr_sub)
    }
    
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
  
  if (verbose) {
    print(df_dimensions_merge)
  }
  
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

#' Get percentage of overlaps for two GRs
#'
#' Only percentage of BinA will be reported
#'  
#' @param A GenomicRanges A
#' @param B GenomicRanges B
#' @param verbose TRUE/FALSE for verbosity
#'
#' @export
get_perc_overlap <- function(A, B, verbose=FALSE) {
  
  library(regioneR)
  
  total_sum <- sum(width(B))
  overlapped <- regioneR::overlapRegions(A=A, B=B) 
  overlapped$type %>% get_vect()
  overlapped['lengthA'] <- overlapped$endA - overlapped$startA
  overlapped['lengthB'] <- overlapped$endB - overlapped$startB
  
  if (verbose) {
    print(overlapped)  
  }
  #sum(subset(overlapped, type=="AinB")$lengthA)/total_sum*100
  sum(subset(overlapped, type=="BinA")$lengthB)/total_sum*100
}

#' Get count of region lengths for each GR
#'  
#' @param A GenomicRanges A
#' @param B GenomicRanges B
#' @param verbose TRUE/FALSE for verbosity
#'
#' @export
stats_GR <- function(given_GR, length2use = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8),
                     names2use = c("1kbp", "10kbp", "100kbp", "1Mbp", "10Mbp", "100Mbp")) {
  
  
  list_summary <- as.list(unclass(summary(given_GR$length)), check.names = FALSE)
  
  for (i in 1:length(length2use)) {
    list_summary[[ names2use[i] ]] <- get_dim.filt(i = given_GR, filt_length=length2use[i])
  }
  
  list_summary
}

#' Get coverage GR or list of GRs
#'  
#' Get coverage for all samples vs. all positions
#'  
#' @param GR_given
#'
#' @export
get_coverage_GR <- function(GR_given) {
  
  library(GenomicRanges)
  library(plyranges)
  
  GRangesList.coord <- GRangesList(GR_given)
  start(GRangesList.coord)
  
  ## get coverage as score
  rle_acro <- GenomicRanges::coverage(GRangesList.coord)
  rle_acro.gr <- as_ranges(rle_acro)
  
  rle_acro.gr
}

#' Get frequency length
#'  
#' @param GR_given
#'
#' @export
freq_length_GR <- function(given_GR) {
  
  length2use = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
  names2use = c("1kbp", "10kbp", "100kbp", "1Mbp", "10Mbp", "100Mbp")
  
  library(DescTools)
  library(reshape2)
  df <- DescTools::Freq(given_GR$length, breaks = length2use)
  df['level'] <- paste0("<", names2use[2:length(names2use)])
  df['sample'] <- given_GR$sample[1]
  
  
  # cum_bp <- c()
  # for (i in names2use[2:length(names2use)]) {
  #  cum_bp <- c(cum_bp, get_sum_bp(given_GR, filt_length = i))
  # }
  
  cum_bp <- c(
    get_sum_bp(given_GR, min_length = 1e3, max_length = 1e4),
    get_sum_bp(given_GR, min_length = 1e4, max_length = 1e5),
    get_sum_bp(given_GR, min_length = 1e5, max_length = 1e6),
    get_sum_bp(given_GR, min_length = 1e6, max_length = 1e7),
    get_sum_bp(given_GR, min_length = 1e7, max_length = 1e8)
  )
  df['bp'] <- cum_bp
  
  df['perc_bp'] <- df$bp/sum(df$bp)
  df['cum_bp'] <- cumsum(df$perc_bp)
  
  return(df)
}

#' Get length of regions for a list of GRs
#'  
#' @param list_GR List of GRs to compute length
#'
#' @export
get_length_list_GRs <- function(list_GR) {
  total_df <- data.frame()
  for (i in names(list_GR)) {
    
    df <- freq_length_GR(list_GR[[i]])
    df$cumfreq <- NULL
    df$cumperc <- NULL
    df$perc <- NULL
    df$perc_bp <- NULL
    df$cum_bp <- NULL
    
    total_df <- rbind(total_df, melt(df))
  }
  total_df$level <- factor(total_df$level, levels=c("<10kbp", "<100kbp", "<1Mbp", "<10Mbp", "<100Mbp"))
  total_df
}



