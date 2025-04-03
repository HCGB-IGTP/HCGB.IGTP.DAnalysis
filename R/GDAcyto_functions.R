get_dim <- function(i) { dim(as.data.frame(i, row.names = NULL))[1] }

# length > 10.000 bp
#' Title
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

#' Title
#'
#' @param matrixinp_given 
#' @param meta_data_given 
#'
#' @export
#'
pheatmap_log <- function(matrixinp_given, meta_data_given=NULL) {
  
  matrixinp.log <- log2(matrixinp_given)
  matrixinp.log[matrixinp.log==-Inf]=0
  
  library(pheatmap)
  pheatmap(matrixinp.log, cluster_rows = TRUE, cluster_cols = TRUE,  
           #scale = "row", 
           annotation_col = meta_data_given)
  
}

#' Get PFB information
#' 
#' Get Population frequency allele-B information for a given population
#'
#' @param variant_data Variant annotation file, provided by illumina
#' @param AF2get Allele frequency to obtain. As default: Global.Allele Freq. from 1k genomes
#'
#' @export
get_pfb <- function(variant_data, AF2get="Global.Allele.Frequency.from.1000G") {
  
  ## VariantAnnotation file avilable in: Illumina oficial website as annotation support files: GDACyto_20047166_A2_VariantAnnotationFile_2022-01-18.tsv
  ## 
  ## /imppc/labs/lslab/share/data/proc_data/20240701_GDACyto_SNP-CNV_acromegalia_JGil_MPuig/data/Illumina_annotation_files/annotation_data/support_files/GDACyto_20047166_A2_VariantAnnotationFile_2022-01-18.tsv
  ## VariantAnnotation_info <- read.csv(VariantAnnotation_info.file, header = TRUE, sep = "\t")

  print("Subsetting:")
  print(AF2get)
  
  snpInfo <- variant_data[,c("Name", "Chr","Position", AF2get)]
  colnames(snpInfo)[4] <-"PFB"
  snpInfo <- na.omit(snpInfo)
  
  print("head(snpInfo)")
  print(head(snpInfo))
  print("dim(snpInfo)")
  print(dim(snpInfo))
  
  # sort by chr position
  library(dplyr)
  snpInfo <- snpInfo %>% arrange(Position)
  
  # set rownames
  row.names(snpInfo) <- snpInfo$Name
  return(snpInfo)
}


#' Extended genoCN plot for BAF and LRR
#'
#' @param pos 
#' @param LRR 
#' @param BAF 
#' @param chr2plot 
#' @param sampleIDs 
#' @param plot_names 
#' @param fileNames 
#' @param types 
#' @param CNA 
#' @param main 
#' @param LRR.ylim 
#' @param cex 
#' @param plot.lowess 
#'
#' @export
#'
plotCN.here <- function (pos, LRR, BAF, chr2plot = NULL, sampleIDs = NULL, 
                         plot_names="genoCN", ## new
                         fileNames = NULL, 
                         types = "genoCN", CNA = TRUE, main = "", LRR.ylim = NULL, 
                         cex = 0.5, plot.lowess = TRUE) {
  
  ## as original from genoCN package
  ## Small modification in the name of each plot provided to easily identify new samples, etc
  ## 
  
  library(genoCN)
  
  if (!is.null(fileNames)) {
    if (length(types) != length(fileNames)) {
      stop("fileNames and types have different lengths\n")
    }
  }
  
  if (length(sampleIDs) == 1) {
    sampleIDs = rep(sampleIDs, length(fileNames))
  }
  if (any(!types %in% c("genoCN", "pennCNV"))) {
    stop("unkown types\n")
  }
  if (any(diff(pos) < 0)) {
    stop("SNPs must be ordered by their positinos\n")
  }
  
  wNA = which(is.na(pos))
  if (!is.null(LRR)) {
    wNA = union(wNA, which(is.na(LRR)))
  }
  if (!is.null(BAF)) {
    wNA = union(wNA, which(is.na(BAF)))
  }
  nplot = 0
  if (!is.null(LRR)) 
    nplot = nplot + 1
  if (!is.null(BAF)) 
    nplot = nplot + 1
  nplot = nplot + length(fileNames)
  if (nplot > 1) {
    par(mfrow = c(nplot, 1))
  }
  if (!is.null(BAF)) {
    par(mar = c(0, 4, 2, 2))
    plot(pos, BAF, xaxt = "n", bty = "n", cex.lab = 1.2, 
         cex.axis = 1.1, main = main, cex = cex)
    abline(h = seq(0, 1, by = 0.25), lty = 2)
  }
  if (!is.null(LRR)) {
    par(mar = c(2, 4, 1, 2))
    plot(pos, LRR, bty = "n", cex.lab = 1.2, cex.axis = 1.1, 
         cex = cex, ylim = LRR.ylim)
    if (plot.lowess) {
      nna = which(!is.na(LRR))
      lines(lowess(LRR[nna] ~ pos[nna], f = 1/50), lwd = 2, col = "red")
    }
    abline(h = 0, lty = 2)
  }
  cols = c("lightblue", "darkblue", "black", "darkgreen", "orange", 
           "darkred", "red", "purple", "brown")
  
  col.pCNV = function(states, cols) {
    colp = character(length(states))
    colp[which(states == 0)] = cols[3]
    colp[which(states == 1)] = cols[4]
    colp[which(states == 2)] = cols[2]
    colp[which(states == 3)] = cols[5]
    colp[which(states == 4)] = cols[7]
    colp
  }
  
  cn.genoCN = function(states, CNA) {
    cn1 = states
    cn1[which(states == 1 | states == 2)] = 2
    cn1[which(states == 3)] = 0
    cn1[which(states == 4)] = 1
    if (CNA) {
      cn1[which(states %in% 5:6)] = 3
      cn1[which(states %in% 7:9)] = 4
    }
    else {
      cn1[which(states == 5)] = 3
      cn1[which(states == 6)] = 4
    }
    cn1
  }
  
  if (length(fileNames) > 0) {
    for (k in 1:length(fileNames)) {
      fileName = fileNames[k]
      type = types[k]
      plot_name=plot_names[k] ## new
      sampleID = sampleIDs[k]
      extSize = file.info(fileName)$size
      par(mar = c(0, 4, 1, 2))
      plot(range(pos[!is.na(pos)]), c(-0.5, 4.5), xaxt = "n", 
           yaxt = "n", bty = "n", type = "n", cex.lab = 1.2, 
           cex.axis = 1.2, xlab = "", ylab = "")
      mtext(0:4, side = 2, line = 0.5, at = 0:4, cex = 0.75)
      mtext(plot_name, side = 2, line = 3, cex = 0.8) ## new
      if (extSize == 0) {
        next
      }
      if (type == "genoCN") {
        ext = read.table(fileName, stringsAsFactors = FALSE, 
                         header = TRUE)
      }
      if (type == "pennCNV") {
        ext = read.table(fileName, stringsAsFactors = FALSE)
        names(ext) = c("chr", "start", "end", "state", "sample", "snp1", "snp2", "score")
      }
      if (!is.null(sampleID)) {
        ext = ext[ext$sample == sampleID, ]
      }
      if (nrow(ext) == 0) {
        start = pos[1]
        end = pos[length(pos)]
        cn1 = 2
        col1 = cols[1]
        rect(start, cn1 - 0.3, end, cn1 + 0.3, density = NULL, 
             angle = 45, col = col1, border = col1, lwd = 0.5)
        next
      }
      chrs = unique(ext$chr)
      if (length(chrs) > 1 & is.null(chr2plot)) {
        stop("there are multiple chromosomes in input file, need to specify which\n         chromosome to plot\n")
      }
      if (!is.null(chr2plot)) {
        ext = ext[ext$chr == chr2plot, ]
      }
      if (nrow(ext) == 0) {
        start = pos[1]
        end = pos[length(pos)]
        cn1 = 2
        col1 = cols[1]
        rect(start, cn1 - 0.3, end, cn1 + 0.3, density = NULL, 
             angle = 45, col = col1, border = col1, lwd = 0.5)
        next
      }
      ext = ext[order(ext$start), ]
      ww = which(pos < ext$start[1])
      if (length(ww) > 0) {
        start = pos[ww[1]]
        end = pos[ww[length(ww)]]
        cn1 = 2
        col1 = cols[1]
        rect(start, cn1 - 0.3, end, cn1 + 0.3, density = NULL, 
             angle = 45, col = col1, border = col1, lwd = 0.5)
      }
      for (i in 1:(nrow(ext) - 1)) {
        start = ext$start[i]
        end = ext$end[i]
        stat = ext$state[i]
        if (type == "genoCN") {
          cn1 = cn.genoCN(stat, CNA)
          col1 = cols[stat]
        }
        if (type == "pennCNV") {
          cn1 = stat
          col1 = col.pCNV(stat, cols)
        }
        rect(start, cn1 - 0.3, end, cn1 + 0.3, density = NULL, 
             angle = 45, col = col1, border = col1, lwd = 0.5)
        ww = which(pos > end & pos < ext$start[i + 1])
        if (length(ww) > 0) {
          start = pos[ww[1]]
          end = pos[ww[length(ww)]]
          cn1 = 2
          col1 = cols[1]
          rect(start, cn1 - 0.3, end, cn1 + 0.3, density = NULL, 
               angle = 45, col = col1, border = col1, lwd = 0.5)
        }
      }
      i = nrow(ext)
      start = ext$start[i]
      end = ext$end[i]
      stat = ext$state[i]
      if (type == "genoCN") {
        cn1 = cn.genoCN(stat, CNA)
        col1 = cols[stat]
      }
      if (type == "pennCNV") {
        cn1 = stat
        col1 = col.pCNV(stat, cols)
      }
      rect(start, cn1 - 0.3, end, cn1 + 0.3, density = NULL, 
           angle = 45, col = col1, border = col1, lwd = 0.5)
      ww = which(pos > end)
      if (length(ww) > 0) {
        start = pos[ww[1]]
        end = pos[ww[length(ww)]]
        cn1 = 2
        col1 = cols[1]
        rect(start, cn1 - 0.3, end, cn1 + 0.3, density = NULL, 
             angle = 45, col = col1, border = col1, lwd = 0.5)
      }
    }
  }
}

#' Extended genoCN plot for BAF and LRR using a Zoom region
#'
#' @param segment_file 
#' @param snpInfo.subset.df 
#' @param sample2process_data 
#' @param start_bp 
#' @param stop_bp 
#' @param several_tracks 
#' @param types.given 
#' @param main_title 
#' @param plot_names.here 
#' @param verbose 
#'
#' @export
#'
plotCN_zoom <- function(segment_file, snpInfo.subset.df, sample2process_data,
                        start_bp=0, stop_bp=NULL, several_tracks=NULL, types.given = c("genoCN","genoCN"),
                        main_title="example plot", plot_names.here=NULL, verbose=FALSE) {
  
  
  ##
  segment_Data <- read.table(segment_file, header = TRUE)
  if (verbose) {
    print("## segment_file data")
    print(segment_file)
    print(head(segment_Data))
  }
  
  if (is.null(stop_bp)) {
    stop_bp <- max(segment_Data$start)
    if (verbose) {
      print("## set stop_bp to max")
    }
  }
  
  
  ## Segment data subset
  seg.subset2check <- subset(segment_Data, start < stop_bp & start > start_bp)
  
  if (verbose) {
    print("## head(seg.subset2check)")
    print(head(seg.subset2check))
  }
  
  snpInfo.subset2check <- subset(snpInfo.subset.df, Position < stop_bp & Position > start_bp)
  sample2process_data2check <- sample2process_data[snpInfo.subset2check$Name,]
  
  if (verbose) {
    print("## head(snpInfo.subset2check)")
    print(head(snpInfo.subset2check))
    print("## head(sample2process_data2check)")
    print(head(sample2process_data2check))
  }
  
  library(genoCN)
  
  if (is.null(several_tracks)) {
    p1 <- plotCN.here(pos=snpInfo.subset2check$Position, LRR=sample2process_data2check$LRR, 
                      BAF=sample2process_data2check$BAF, 
                      main = main_title, fileNames=segment_file)
    
  } else {
    if (verbose) {
      print("## Several tracks provided")
      print(c(segment_file, several_tracks))
    }      
    
    p1 <- plotCN.here(pos=snpInfo.subset2check$Position, LRR=sample2process_data2check$LRR, 
                      BAF=sample2process_data2check$BAF, types=types.given, plot_names = plot_names.here,
                      main = main_title, fileNames=c(segment_file, several_tracks))
  }
  
  print(p1)
  
  plotCN_zoom.list <- list(
    "data" = seg.subset2check,
    "plot" = p1
  ) 
  
  return(plotCN_zoom.list) 
}


#' Create GR list subseting dataframe
#'
#' @param data2use 
#'
#' @return
#' @export
#'
#' @examples
GR_list_samples <- function(data2use) {
  GR_list <- list()
  for (sample_name in levels(as.factor(data2use[['sample']]))) {
    print(paste0("+ Subsetting: ", sample_name, "...."))
    tmp.gr <- GenomicRanges::makeGRangesFromDataFrame(df = subset(data2use,sample==sample_name),  
                                                      keep.extra.columns = TRUE, na.rm = TRUE)
    
    GR_list[[sample_name]] = tmp.gr
  } 
  print(" + Done!")
  print("")
  return(GR_list)
}


#' Read pennCNV output file
#'
#' @param PennCNV_file 
#'
#' @return
#' @export
#'
#' @examples
read_pennCNV <- function(PennCNV_file) {
  
  data_pennCNV <- read.table(file = PennCNV_file)
  names(data_pennCNV) <- c("chr", "start", "end", "state", 
                           "sample", "snp1", "snp2", "score", "num_snp")
  
  ## add new columns
  data_pennCNV['id'] <- paste0(data_pennCNV$chr, "_", 
                               data_pennCNV$start, "_",
                               data_pennCNV$end)
  data_pennCNV['length'] <- data_pennCNV$end - data_pennCNV$start
  data_pennCNV['software'] <- "PennCNV"
  
  print("head(data_pennCNV)")
  print(head(data_pennCNV))
  
  ## return dataframe and GR list for each sample
  list_data_pennCNV = list(
    "df"=data_pennCNV,
    "GR.list" = GR_list_samples(data2use = data_pennCNV)
  )
  
  
  return(list_data_pennCNV)
}

#' Read genoCN output file
#'
#' @param genoCN_file 
#' @export
read_genoCN <- function(genoCN_file) {
  
  genoCN_data <- read.table(file = genoCN_file, header = TRUE)
  genoCN_data['id'] <- paste0(genoCN_data$chr, "_", genoCN_data$start, "_", genoCN_data$end)
  genoCN_data['length'] <- genoCN_data$end - genoCN_data$start
  genoCN_data['software'] <- "genoCN"
  colnames(genoCN_data)[10] <- 'num_snp'
  
  print("head(genoCN_data)")
  print(head(genoCN_data))
  
  ## return dataframe and GR list for each sample
  list_data_genoCN = list(
    "df"=genoCN_data,
    "GR.list" = GR_list_samples(data2use = genoCN_data)
  )
  
  
  return(list_data_genoCN)
}

