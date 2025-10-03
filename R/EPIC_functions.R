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

#' Get limma results for EPIC probes
#'
#' @param normData.given Dataframe with probes as rownames
#' @param sample_sheet.given Dataframe with information
#' @param design.given Statistical design
#' @param comp.name.given Name for the comparison variable
#' @param folder_given Folder to store results
#' @param contrasts.matrix.here Contrasts matrix to test
#' @param illumina_annot Illumina annotation dataframe
#' @param b_values_given Dataframe containing beta_values for each probe and sample (columns)
#' @param FCcutoff.given Log FC difference (absolute values) for filtering
#' @param mean.diff Mean difference (absolute values) for filtering
#'
#' @export
get_limma_EPIC_results <- function(normData.given, sample_sheet.given, design.given, comp.name.given, 
                                   folder_given, contrasts.matrix.here, illumina_annot, 
                                   b_values_given, FCcutoff.given = log2(1.2), diff.mean.given=0.2, paired_variable=NULL) {
  
  ## create limma comparison
  res_limma <- get_limma_results(normData = normData.given,
                                 df_treatment_Ind = sample_sheet.given,
                                 design.given = design.given,
                                 comp_name = comp.name.given, 
                                 FCcutoff.given = FCcutoff.given,
                                 contrasts.matrix.given = contrasts.matrix.here,
                                 coef.given.list = colnames(contrasts.matrix.here),
                                 OUTPUT_Data_sample = folder_given,
                                 EPIC=TRUE, 
                                 forceResults = TRUE, 
                                 paired_variable=paired_variable)
  
  save(res_limma, file=file.path(folder_given, paste0(comp.name.given, ".RData")))
  
  ## annotate and return information 
  GR_list <- list()
  GR_filt <- list()
  Samplslist <- list()
  for (x in names(res_limma$results_list)) {
    
    print(paste0("+ Comparison", x))
    result_names_list <- get_comparison_resultsNames(str_given = paste0(comp.name.given, "_", x))
    print(result_names_list)                                                 
    
    all_Data.res <- res_limma$results_list[[x]][['alldata.norm.res']]
    rownames(all_Data.res) <- all_Data.res$Gene
    
    ## get beta values
    all_Data.res <- get_beta_mean_vals(table_dat = all_Data.res, beta_val = b_values_given, 
                            list1 = res_limma$results_list[[x]][['Samplslist']][['numerator']],
                            list2 = res_limma$results_list[[x]][['Samplslist']][['denominator']])
    ##
    print("+ Candidate probes")
    print(length(rownames(subset(all_Data.res, diff.mean>diff.mean.given & adj.P.Val<0.05))))
    
    ## merge with annotation 
    all_Data.res <- merge(x = illumina_annot[rownames(all_Data.res),
                                           c("CpGchrm", "probeStrand", "CpGbeg", "CpGend")],
                          y = all_Data.res, by="row.names")
    
    ## create genomic ranges
    all_Data.norm.GR <- GenomicRanges::makeGRangesFromDataFrame(all_Data.res,  
                                          start.field = "CpGbeg", end.field = "CpGend", 
                                          seqnames.field =  "CpGchrm", strand.field = "probeStrand",
                                          keep.extra.columns = TRUE)
    GR_list[[x]] <- all_Data.norm.GR
    
    GR_filt[[x]] <- subset(all_Data.norm.GR, diff.mean>0.2 & adj.P.Val<0.05)
    
    Samplslist[[x]] <- res_limma$results_list[[x]][['Samplslist']]
    
  }
  
  return(list("Samplslist" = Samplslist,
              "GR_list" = GR_list,
              "GR_filt"=GR_filt))
}


#' Call DMRcate for several comparisons
#'
#' @param sample_metadata.given Dataframe with information
#' @param var_name Name for the comparison variable
#' @param dir_folder.given Folder to store results
#' @param m_values.noSNPs.given Dataframe with probes as rownames, m values with no SNPs (as a result of DMRcate::rmSNPandCH)
#' @param lambda.given dmrcate parameter. Gaussian kernel bandwidth for smoothed-function estimation. Default: 1000
#' @param C.given dmrcate parameter. Scaling factor for bandwidth. Gaussian kernel is calculated where lambda/C = sigma. Default: 2
#' @param fdr.given False discovery rate cutoff. Default 0.25.
#'
#' @export
call_dmrcate_comparisons <- function(sample_metadata.given, var_name,  dir_folder.given, m_values.noSNPs.given,
                         lambda.given=1000, C.given=2, fdr.given=0.05, paired_variable=NULL, other_var_name=NULL) {
  
  library(DMRcate)
  library(ExperimentHub)
  library(minfi)
  library(SummarizedExperiment)
  
  dir.create(dir_folder.given)
  
  possible_comp <- get_possible_comparisons(sample_metadata.given[[var_name]])
  original_test.df <- sample_metadata.given
  
  for (i in 1:length(possible_comp)) {
    ref_level <- unlist(possible_comp[i])[1]
    comp_level <- unlist(possible_comp[i])[2]
    
    print("---------------------------------")
    my_cmp <- paste0(comp_level, "_vs_", ref_level)
    print(my_cmp)
    print("---------------------------------")
    
    ##---------------------------------
    ## check if already generated
    ##---------------------------------
    my_dir <- file.path(dir_folder.given, my_cmp)
    dir.create(my_dir)
    
    if (file.exists(file.path(my_dir, "dmrcate_results.RData"))) {
      print("++ Results already generated")
      print("")
      next
    }
    ##---------------------------------
    
    ##---------------------------------
    ## Prepare
    ##---------------------------------
    ## subset only rows of interest
    test.df <- original_test.df[original_test.df[[var_name]] %in% c(ref_level, comp_level),]
    
    ## set reference
    levels_var <- levels(as.factor(test.df[[var_name]]))
    test.df[[var_name]] <- factor(test.df[[var_name]], levels = c(ref_level, comp_level))
    
    print(paste0("+ Set reference level: ", ref_level))
    print(test.df[[var_name]])
    ##---------------------------------
    
    ##---------------------------------
    ## create model matrix
    ##---------------------------------
    print("+ Create model matrix")
    if (is.null(other_var_name)) {
      design.cmp <- model.matrix(~test.df[[var_name]]) # describe model to be fit
    } else {
      design.cmp <- model.matrix(~test.df[[var_name]]+test.df[[other_var_name]]) # describe model to be fit
      colnames(design.cmp)[3] <- other_var_name 
      }
    colnames(design.cmp)[2] <- comp_level  
    
    print(head(design.cmp))
    ##---------------------------------
    
    ##---------------------------------
    ## call DMRcate for each comparison generated
    ##---------------------------------
    call_dmrcate(m_values.noSNPs.given = m_values.noSNPs.given, coef.given=2,
                 metadata.given=test.df, design.cmp = design.cmp,
                 my_cmp=my_cmp, my_dir=my_dir,
                 lambda.given = lambda.given,  
                 C.given = C.given, fdr.given = fdr.given, paired_variable = paired_variable)
    ##---------------------------------
  } 
}


#' Call DMRcate and previous steps
#'
#' @param m_values.noSNPs.given 
#' @param my_cmp 
#' @param my_dir 
#' @param metadata.given 
#' @param design.cmp 
#' @param fdr.given 
#' @param coef.given 
#' @param lambda.given 
#' @param C.given 
#' @param paired_variable
#'
#' @export
call_dmrcate <- function(m_values.noSNPs.given, my_cmp, my_dir, 
                         metadata.given, design.cmp, fdr.given=0.05, coef.given=2, lambda.given=1000, C.given=2, paired_variable=NULL) {
  
  
  library(DMRcate)
  library(ExperimentHub)
  library(minfi)
  library(SummarizedExperiment)
  
  ##---------------------------------
  # CpG annotate
  ##---------------------------------
  print(paste0("+ Call DMRcate for ", my_cmp))
  print("+ CpG annotate: ")
  myannotation.cmp <- cpg.annotate("array", 
                                   object=m_values.noSNPs.given[,rownames(metadata.given)], 
                                   what = "M", arraytype = "EPICv2", epicv2Filter = "mean",
                                   epicv2Remap = TRUE, 
                                   analysis.type="differential",
                                   design=design.cmp, 
                                   fdr = fdr.given, block=paired_variable,
                                   coef=coef.given)
  print("myannotation.cmp")
  print(myannotation.cmp)
  print(class(myannotation.cmp))
  ##---------------------------------
  
  ##---------------------------------
  # dmrcate
  ##---------------------------------
  print("+ Call dmrcate: ")
  dmrcoutput.cmp <- try({
    dmrcate(myannotation.cmp, lambda=lambda.given, C=C.given)
  }, silent = TRUE)
  
  # Process any error messages
  if (!is(dmrcoutput.cmp, "DMResults")) {
    
    print("+ CpG annotate (less restrictive): ")
    myannotation.cmp <- cpg.annotate("array", 
                                     object=m_values.noSNPs.given[,rownames(metadata.given)], 
                                     what = "M",
                                     arraytype = "EPICv2", 
                                     epicv2Filter = "mean",
                                     epicv2Remap = TRUE, 
                                     analysis.type="differential",
                                     design=design.cmp, 
                                     fdr = 1, 
                                     coef=coef.given)
    
    dmrcoutput.cmp <- try({
      dmrcate(myannotation.cmp, lambda=lambda.given, C=C.given)
    }, silent = TRUE)
    
    # Process any error messages
    if (!is(dmrcoutput.cmp, "DMResults")) {
      print("+ No DMRs regions found for this comparison")
      save(dmrcoutput.cmp, myannotation.cmp,
           file = file.path(my_dir, "dmrcate_results.RData"))
      return()
    }
  }
  
  print("dmrcoutput.cmp")
  print(dmrcoutput.cmp)
  print(class(dmrcoutput.cmp))
  ##---------------------------------
  
  ##---------------------------------
  ## Get summary of results
  ##---------------------------------
  print("+ Get ranges: ")
  results.ranges.cmp <- try({ 
    extractRanges(dmrcoutput.cmp, genome = "hg38")}, 
    silent = TRUE)
  
  if (!is(results.ranges.cmp, "GRanges")) {
    print("+ No DMRs regions found for this comparison")
    save(dmrcoutput.cmp, myannotation.cmp,
         file = file.path(my_dir, "dmrcate_results.RData"))
    return()
  } 
  
  results.ranges.cmp$abs_meandiff <- abs(results.ranges.cmp$meandiff)
  results.ranges.cmp <- results.ranges.cmp[order(results.ranges.cmp$abs_meandiff, decreasing = TRUE),]
  
  print("head(results.ranges)") 
  print(head(results.ranges.cmp)) 
  
  print("+ Subset significant results") 
  results.ranges.filt <- subset(results.ranges.cmp, min_smoothed_fdr<0.05 & abs_meandiff>0.2)
  ##---------------------------------
  
  ##---------------------------------
  ## Save
  ##---------------------------------
  print("++ Save results!")
  try(write.csv(x = GenomicRanges::as.data.frame(results.ranges.cmp),  file = file.path(my_dir, "dmrcate_results.all.csv")))
  try(write.csv(x = GenomicRanges::as.data.frame(results.ranges.filt), file = file.path(my_dir, "dmrcate_results.sign.csv")))
  
  save(results.ranges.filt, results.ranges.cmp, dmrcoutput.cmp, myannotation.cmp, 
       file = file.path(my_dir, "dmrcate_results.RData"))
}

#' Call limma for several comparisons
#'
#' @param sample_metadata.given Dataframe with information
#' @param var_name Name for the comparison variable
#' @param dir_folder.given Folder to store results
#' @param values_df Dataframe with probes as rownames, m values with no SNPs (as a result of DMRcate::rmSNPandCH)
#' @param fdr.given False discovery rate cutoff. Default 0.05.
#' @param other_var_name Other variable to include in the design
#' @param paired_variable If pair samples, include this variable
#' @param EPICannotation.given EPIC annotation dataframe
#' @param bvalues Beta values dataframe 
#'
#' @export
call_limma_comparisons <- function(sample_metadata.given, var_name,  dir_folder.given, values_df,
                                   fdr.given=0.05, EPICannotation.given, bvalues,
                                   paired_variable=NULL, other_var_name=NULL) {
  
  
  dir.create(dir_folder.given)
  
  possible_comp <- get_possible_comparisons(sample_metadata.given[[var_name]])
  original_test.df <- sample_metadata.given
  
  for (i in 1:length(possible_comp)) {
    ref_level <- unlist(possible_comp[i])[1]
    comp_level <- unlist(possible_comp[i])[2]
    
    print("---------------------------------")
    my_cmp <- paste0(comp_level, "_vs_", ref_level)
    print(my_cmp)
    print("---------------------------------")
    
    ##---------------------------------
    ## check if already generated
    ##---------------------------------
    my_dir <- file.path(dir_folder.given, my_cmp)
    dir.create(my_dir)
    
    if (file.exists(file.path(my_dir, "data2return.RData"))) {
      print("++ Results already generated")
      print("")
      next
    }
    ##---------------------------------
    
    ##---------------------------------
    ## Prepare
    ##---------------------------------
    ## subset only rows of interest
    test.df <- original_test.df[original_test.df[[var_name]] %in% c(ref_level, comp_level),]
    
    ## set reference
    levels_var <- levels(as.factor(test.df[[var_name]]))
    test.df[[var_name]] <- factor(test.df[[var_name]], levels = c(ref_level, comp_level))
    
    print(paste0("+ Set reference level: ", ref_level))
    print(test.df[[var_name]])
    ##---------------------------------
    
    ##---------------------------------
    ## create model matrix
    ##---------------------------------
    print("+ Create model matrix")
    if (is.null(other_var_name)) {
      design.cmp <- model.matrix(~test.df[[var_name]]) # describe model to be fit
    } else {
      design.cmp <- model.matrix(~test.df[[var_name]]+test.df[[other_var_name]]) # describe model to be fit
      colnames(design.cmp)[3] <- other_var_name 
    }
    colnames(design.cmp)[2] <- comp_level  
    
    print(head(design.cmp))
    
    ## create contrast
    contrasts.matrix.cmp <- make_all_contrasts(group = c(comp_level, ref_level))
    ##---------------------------------
    
    ##---------------------------------
    ## call limma for each comparison generated
    ##---------------------------------
    
    res.limma <- try(get_limma_EPIC_results(normData.given = values_df[,rownames(test.df)]),
        sample_sheet.given = test.df,
        design.given = design.cmp,
        comp.name.given = var_name, 
        contrasts.matrix.here = contrasts.matrix.cmp,
        folder_given = my_dir, 
        illumina_annot = EPICannotation.given, b_values_given = bvalues)

    if (is.null(res.limma)) {
      
    } else {
      save(res.limma, file=file.path(my_dir, paste0(var_name, ".RData")))
    }

  } 
}

