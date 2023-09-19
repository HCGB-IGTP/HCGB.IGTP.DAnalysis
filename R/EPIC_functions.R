#' Annotates GR Report for each region of bumps provided
#'
#' This functions create a GR object with the regions of interest from bump_table provided
#' @param bump_table Bumphunter dataframe
#' @param out_report Not used anymore, to be discarded
#' @param df_chr Dataframe generated as in get_chr_len() containing names and length for each Chromosome of interest
#' @param EPIC Set v1 or v2 for annotation. It uses hg19 for EPICv1 or hg38 for EPICv2. Default: v2
#' @export
create_report <- function(bump_table, out_report, df_chr, EPIC='v2') {
  library(GenomicRanges)
  
  ## Use hg19 for EPICv1 or hg38 for EPICv2   
  if (EPIC=="v2") {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    
  } else {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  
  ####
  regions <- GRanges(
    seqnames = bump_table$chr,
    IRanges(start = bump_table$start, end=bump_table$end),
    strand="*", 
    value=bump_table$value, 
    area=bump_table$area,
    cluster=bump_table$cluster, 
    L=bump_table$L, 
    clusterL=bump_table$clusterL,
    
    ## add new
    p.value=bump_table$p.value,
    fwer=bump_table$fwer,
    p.valueArea=bump_table$p.valueArea,
    fwerArea=bump_table$fwerArea
  )
  
  seqlengths(regions) <- df_chr[unique(bump_table$chr),]
  
  #library(regionReport)
  #report <- renderReport(regions, name, pvalueVars=NULL, 
  #densityVars=c("Area" = "area", "Value"="value","Cluster Length" = "clusterL"),significantVar=NULL,output="report", outdir=out_report,device="png"
  #)
  
  genes <- annotateTranscripts(txdb=txdb)  
  annotation <- matchGenes(x = regions, subject = genes)
  if (nrow(annotation) == 0) {
    regions.df <- as.data.frame(regions)  
  } else {
    ## Add annotation information
    regions.df <- cbind(as.data.frame(regions), annotation)
  }
  
  list2return <- list(
    "regions" = regions,
    "regions.df" = regions.df
  )
  
  return(list2return)
}


#' Filter bumps and retain more significant
#'
#' This functions filters a bump_table provided. It uses only the cluster L given. It could also use fwer as a method to discard false positives
#' @param bump_table Bumphunter dataframe
#' @param clusterL_given Number of probes to be within each cluster
#' @param L_given Numberof probes with the same trending line within each cluster
#' @export
filter_bump_table <- function(bump_table, clusterL_given=3, L_given=3) {
  
  ## create percentage
  bump_table$percL <- (bump_table$L/bump_table$clusterL)*100
  
  ## filter
  bump_table2 <- bump_table[bump_table$clusterL>=clusterL_given,]
  bump_table3 <- bump_table2[bump_table2$L>=L_given,]
  
  ## sort
  bump_table3 <- bump_table3[order(bump_table3$percL, decreasing = TRUE),]
  
  return(bump_table3)
}

#' Save BUMP information in BED file format
#'
#' This functions saves a bump_table provided in BED format for later analysis
#' @param bump_table Bumphunter dataframe
#' @param file_out File to store results
#' @export
save_bed_info <- function(bump_table, file_out) {
  
  ## probe annotation bed information
  bump_table$ID <- paste(bump_table$chr, bump_table$start, sep="_")
  bump_table_sort <- bump_table[order(bump_table$chr),]
  write.table(x=bump_table_sort[,c(1,2,3,4,5,9,10,11,12,15)], file=file_out, 
              sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}


#' Return dataframe with chromosome length for hg19 genome
#'
#' This functions returns chromosome length for hg19 genome
#' @param genome_version Set hg38 or hg19. Default: hg38
#' @export
get_chr_len <- function(genome_version='hg38'){
  
  ## hg38 information
  if (genome_version=='hg38') {
    library(BSgenome.Hsapiens.UCSC.hg38)
    print("hg38 information")
    df_chr <- as.data.frame(seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))[1:25,]
    df_chr$isCircular <- NULL
    df_chr$genome <- NULL
    colnames(df_chr)[1] <- 'chr_length'
    
    ## hg19 information
  } else if (genome_version=='hg19') {
    print("hg19 information")
    ## get chr length
    chr_list <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                  "chr20","chr21","chr22")
    chr_length <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,
                    135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,
                    59128983,63025520,48129895,51304566)
    
    df_chr <- data.frame(chr_list, chr_length)
    rownames(df_chr) <- df_chr$chr_list
    df_chr$chr_list <- NULL
  } else {
    print("Not implemented yet")
    df_chr <- data.frame()
  }
  
  return(df_chr)
}

#' Create limma analysis for each probe
#'
#' This functions creates a simplified version of limma analysis for each probe
#' @param values_df Values to test
#' @param design_given Statistical design
#' @param contrasts.matrix_given Contrasts to test
#' @param pval_sign Cufoff for pvalue adjusted
#' @param thr Threshold for LFC
#' @export
get_probes_limma <- function(values_df, design_given, contrasts.matrix_given , pval_sign=0.05, thr=0.26) {
  ## linear model regression fit analysis
  print("Linear model regression fit analysis")
  print("Fit each probeset to model")
  fit_1 <- limma::lmFit(values_df, design_given)  # fit each probeset to model
  
  print("Create contrast fit")
  fit_2 <- limma::contrasts.fit(fit_1, contrasts.matrix_given)
  
  print("Empirical Bayes adjustment")
  efit_3 <- limma::eBayes(fit_2, trend=TRUE)        # empirical Bayes adjustment
  
  print("Decide Tests")
  results_efit <- decideTests(efit_3)
  
  print("TopTable")
  res_limma <- limma::topTable(efit_3, number= 5e6, adjust.method = "BH")
  
  print("Delete NA")
  res_limma <- delete_na(res_limma)
  res_limma <- res_limma[order(res_limma$B, decreasing = TRUE),]
  
  ## Filter significant
  res_limma_filter <- res_limma[res_limma$adj.P.Val<pval_sign & abs(res_limma$logFC)>thr,]
  
  results_limma <- list(
    'logFC' = thr,
    'adj.P.Value'= pval_sign,
    'res_limma_filter' = res_limma_filter,
    'res_limma' = res_limma,
    'fit_object' = fit_2,
    'efit_object' = efit_3
  )
  
  return(results_limma)
}

#' Get mean beta values for two list of samples
#'
#' This functions obtains a mean beta value for each group of samples
#' @param table_dat Dataframe with probes as rownames
#' @param beta_val Dataframe containing beta_values for each probe and sample (columns)
#' @param list1 Vector of samples to test fro group1
#' @param list1 Vector of samples to test fro group2
#' @export
get_beta_mean_vals <- function(table_dat, beta_val, list1, list2) {
  ## Calculate diff mean
  table_dat['mean.group1'] <-rowMeans(beta_val[rownames(table_dat),list1])
  table_dat['mean.group2'] <-rowMeans(beta_val[rownames(table_dat),list2])
  table_dat['diff.mean'] <- abs(table_dat$mean.group1 - table_dat$mean.group2)
  table_dat['sign'] <- ifelse(table_dat$mean.group1 - table_dat$mean.group2 > 0,"up", "down")
  
  return(table_dat)
}


#' Filter dmpFinder methylation stats
#'
#' This functions obtains a mean beta value for each group of samples and filters based on statistics
#' @param dmpData Dataframe with probes as rownames
#' @param beta_val Dataframe containing beta_values for each probe and sample (columns)
#' @param list1 Vector of samples to test fro group1
#' @param list1 Vector of samples to test fro group2
#' @param pval Pvalue filtering cutoff. Default: 0.05
#' @param qval Pvalue adjusted filtering cutoff. Default: 0.05
#' @param diff.mean Mean difference (absolute values) for filtering
#' @export
filter_dmpFinder_table <- function(dmpData, beta_val, list1, list2, pval=0.05, qval=0.05, diff.mean=0.2) {
  
  ## get info
  #dmpData <- read.csv(dmpFile, header = TRUE, dec = ".", row.names = 1)
  #dmpData <- na.omit(dmpData)
  
  ## Calculate diff mean
  dmpData <- get_beta_mean_vals(table_dat = dmpData, beta_val = beta_val, list1=list1, list2=list2)
  
  ## filter
  dmpData <- dmpData[dmpData$pval <= pval,]
  dmpData <- dmpData[dmpData$qval <= qval,]
  dmpData <- dmpData[ abs(dmpData$diff.mean) >= diff.mean,]
  
  ## sort
  dmpData <- dmpData[order(dmpData$diff.mean, decreasing = TRUE),]
  
  return(dmpData)
}

