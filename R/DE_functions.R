#' Create MA plot
#'
#' This functions creates MA plots given a dataframe containing:
#' baseMean, log2FoldChange and padj header names.
#' @param res dataframe 
#' @param thresh pvalue adjusted threshold 
#' @param textcex determines text size
#' @export
maplot <- function (res, thresh=0.05, textcex=1, ...) {
  ## Susana MAPLOT
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
}

#'
#' This functions discard low counts in a DESeq object.
#' Uses 10 as default.
#' 
#' @param dds_object DESeq2 object 
#' @param min_count integer for minimum count cutoff [Default=10]. 
#' @export
##  #####
discard_lowCounts = function(dds_object, min_count=10) {
  keep <- rowSums(counts(dds_object)) >= min_count
  return(dds_object[keep,])
}

#' Create Volcano plot
#'
#' This functions creates volcano plots given a dataframe containing:
#' baseMean, log2FoldChange and padj header names.
#' @param res dataframe 
#' @param lfcthresh Log fold change threshold 
#' @param sigthresh pvalue adjusted threshold 
#' @param main Main plot legend 
#' @param legendpos Legend position
#' @export
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|FC|>",round(2^lfcthresh,1),sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#' Discard empty counts
#'
#' This functions prepares countsfile for DESeq discarding rows with
#' many 0 for all samples. It uses the given cutoff, to discard rows containing
#' greater than a cutoff in percentage of 0's. Default cutoff value is 0.9
#' @param countsF Count file for DESeq
#' @param cutoff Cutoff for discarding rows
#' @export
discard_0_counts <- function(countsF, cutoff=0.9) {
  row.0 <- rowSums(countsF==0) # for each gene, number of samples with 0 counts
  row.0_left <- row.0[row.0 <= cutoff * ncol(countsF) ]
  countsFile2 <-countsF[names(row.0_left),]
  return(countsFile2)
}

#' Plot DESeq2 p-values
#'
#' This functions plots pvalues original and after filtering
#' @param pdf_name File path for the pdf to save images
#' @param res Original results (dds object)
#' @param res_filtered Filtered results (dds object)
#' @export
plot_DESeq2_pvalues <- function(pdf_name, res, res_filtered) {
  pdf(pdf_name)
  p1 <- hist(res$pvalue, breaks=100, col="grey50", main = "p-value distribution (No filtered)")
  p2 <- hist(res_filtered$pvalue, breaks=100, col="grey50", main = "p-value distribution (Filtered)")
  
  ### p-value histogram
  use <- res_filtered$baseMean > metadata(res_filtered)$filterThreshold
  h1 <- hist(res_filtered$pvalue[!use], breaks=100, plot=FALSE)
  h2 <- hist(res_filtered$pvalue[use], breaks=100, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  
  plot(p1)
  plot(p2, add=TRUE, col="powderblue")
  plot(h1, add=TRUE, col="khaki")
  legend("top", legend=c("Filtered low count", "Filtered FDR", "Passed"), fill = c("grey", "khaki", "powderblue"))
  dev.off()
}

#' Discard low counts in dataframe
#'
#' This functions discard low counts in a dataframe with any non-numeric column
#' @param df_given Dataframe to pase
#' @param min_count Minimum count to use as cutoff [Default: 10].
#' @export
discard_lowCounts_df = function(df_given, min_count=10) {
  df_given <- df_given %>% mutate(total=rowSums(select_if(., is.numeric)))
  keep <- df_given['total'] >= min_count
  df_given['total'] <- NULL
  return(df_given[keep,])
}

#' Plot DESeq2 p-values
#'
#' This functions plots pvalues original and after filtering
#' @param dds_object DESeq2 object
#' @param coef_n Coefficient number obtain using resultsName(dds)
#' @param name Name of the comparison
#' @param numerator Name of the numerator comparison
#' @param denominator Name of the denominator comparison
#' @param OUTPUT_Data_dir Folder path to store results
#' @param df_treatment_Ind Dataframe containing additional information for each sample
#' @param threads Number of CPUs to use [Default: 2].
#' @export

DESeq2_HCGB_function = function(dds_object, coef_n, comp_name, comp_ID="comp1",
                                numerator="example1", denominator="example2", 
                                OUTPUT_Data_dir, df_treatment_Ind, threads=2, forceResults=FALSE) {
  
  library(DESeq2)
  library(vsn)
  library(EnhancedVolcano)
  library(BiocParallel)
  library(ggfortify)
  library(ggrepel)
  library(pheatmap)
  library(reshape2)
  library(RColorBrewer)
  
  
  ## set name
  file_name <- paste0(comp_ID, "_", comp_name, "_", numerator, "_vs_", denominator)
  
  ## start
  print (paste0("## Starting: ", file_name))
  
  ## Set parallel threads
  print (paste0("Set Multicore: ", as.numeric(threads)))
  register(MulticoreParam(as.numeric(threads)))
  
  ## Create folder
  OUTPUT_Data_sample = file.path(OUTPUT_Data_dir, file_name)
  print (paste0("Create folder: ", OUTPUT_Data_sample))
  dir.create(OUTPUT_Data_sample)
  
  if (forceResults) {
    
  } else {
    ## check if previously done
    if (file.exists(file.path(OUTPUT_Data_sample, "data2return.RData"))) {
      print("Data already available in:")
      print(OUTPUT_Data_sample)
      return()
    }
  }
  
  
  ######################################################################
  ## Generate results according to comparison
  ######################################################################
  
  ## Compare pvalues distribution
  res <- DESeq2::lfcShrink(dds_object, coef=coef_n, type="apeglm", parallel = TRUE)
  
  ## filter low counts
  dds_object = HCGB.IGTP.DAnalysis::discard_lowCounts(dds_object = dds_object)
  res_filtered <- DESeq2::lfcShrink(dds_object, coef=coef_n, type="apeglm", parallel = TRUE)
  
  ######################################################################  
  ### Plot p-values
  ######################################################################
  ## original/filtered
  Pvalues_pdf <- file.path(OUTPUT_Data_sample, paste0(comp_name,"_", numerator, "_vs_", denominator," p-value_distribution.pdf"))
  HCGB.IGTP.DAnalysis::plot_DESeq2_pvalues(Pvalues_pdf, res, res_filtered)
  
  ######################################################################
  # Write Results tables
  ######################################################################
  ## write.table(res_filtered, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_table.txt")), sep="\t", row.names=T, col.names=NA, quote=F)
  
  # Results table ordered by adjusted p-value
  resOrdered <- res_filtered[order(res_filtered$padj),]
  ## write.table(resOrdered, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_padj-ordered_table.txt")), sep="\t", row.names=T, col.names=NA, quote=F)
  
  # Save normalized values
  # Normalized values
  normValues <- counts(dds_object, normalized=T)
  ## write.table(normValues, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_NormValues_table.txt")), sep="\t", row.names=T, col.names=T, quote=F)
  
  ## Merge normalized values and differential expression
  alldata <- merge(as.data.frame(counts(dds_object, normalized=TRUE)), as.data.frame(res_filtered), by="row.names", sort=FALSE)
  names(alldata)[1] <- "Gene"
  ## Write results
  write.table(alldata, file=file.path(OUTPUT_Data_sample, 
                                      paste0(file_name, "-ResultsCounting_NormValues-and-DE_table.txt")), sep="\t", 
              row.names=T, col.names=NA, quote=T)
  
  ##get significant data only with counts in all samples
  ##get significant data only padj < 0.05 and log2FoldChange >1.2
  #sign.data<-alldata[alldata$padj<0.05,]
  #sign.data<-sign.data[!is.na(sign.data$padj),]
  #sign.data<-sign.data[abs(  sign.data$log2FoldChange)>log2(1.2),]
  
  sign.data <- filter_signficant_DESEQ(alldata, sign_value = 0.05, LFC = log2(1.2))
  row.names(sign.data) <- sign.data$Gene
  
  # Results table in the same order than counting table
  write.table(sign.data, 
              file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_table_SignificantDE.txt")), 
              sep="\t", row.names=T, col.names=NA, quote=TRUE)
  
  ## check if it is worth to continue, avoid error if missing sign.data
  print(length(rownames(sign.data)))
  
  if (length(rownames(sign.data)) < 2) {
    
    data2save<- list(
      "alldata" = alldata,
      #"dds_object" = dds_object,
      #"res"=res,
      "res_filtered" = res_filtered,
      "sign.df"=sign.data,
      "sign.genes"=sign.data$Gene,
      "sign.count"=length(sign.data$Gene)
    )
    
    ## dump in disk RData
    save(data2save, file=file.path(OUTPUT_Data_sample, "data2return.RData"))
    
    data2return<- list(
      "res_filtered" = res_filtered,
      "sign.df"=sign.data,
      "sign.genes"=sign.data$Gene,
      "sign.count"=length(sign.data$Gene)
    )
    
  } else {
    
    ######################################################################
    ## ma plot
    ######################################################################
    pdf(file.path(OUTPUT_Data_sample, "DiffExpression-maplot.jpeg"))
    plotMA(res_filtered)
    dev.off()
    
    ######################################################################
    ## volcano
    ######################################################################
    #jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-volcano-plot.jpeg")), 1500, 1000, pointsize=20)
    volcano_main_title <- paste0("Volcano Plot: ", numerator, " vs ", denominator)
    volcan_plot <- EnhancedVolcano::EnhancedVolcano(alldata, x="log2FoldChange", y="padj", lab="",
                                                    pCutoff=0.05, FCcutoff=log2(1.2), pointSize=3, labSize=6) + 
      ggplot2::scale_x_continuous() + ggplot2::labs(title = volcano_main_title)
    
    HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                  name_file = paste0(file_name, "_DiffExpression-volcano-plot"), plot_given = volcan_plot)
    
    ######################################################################
    ## Transform normal
    ######################################################################
    
    ## rlog transformation
    rld <- DESeq2::rlogTransformation(dds_object)
    
    # variance stabilizing
    vsd <- DESeq2::varianceStabilizingTransformation(dds_object, blind = FALSE)
    
    ## The figure below plots the standard deviation of the transformed data, 
    ## across samples, against the mean, using the shifted logarithm transformation (ntd), 
    ## the regularized log transformation (rld) and the variance stabilizing transformation (vst).
    ## The shifted logarithm has elevated standard deviation in the lower count 
    ## range, and the regularized log to a lesser extent, while for the variance 
    ## stabilized data the standard deviation is roughly constant along the whole
    ## dynamic range.
    ## 
    ## Note that the vertical axis in such plots is the square root of the variance 
    ## over all samples, so including the variance due to the experimental conditions. 
    ## While a flat curve of the square root of variance over the mean may seem like 
    ## the goal of such transformations, this may be unreasonable in the case of 
    ## datasets with many true differences due to the experimental conditions.
    
    ## ntd <- normTransform(dds_object)
    ## pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-Transformation.pdf")))
    ## Normal transformation  
    ## meanSdPlot(assay(ntd))
    ## RLE transformation
    ## meanSdPlot(assay(rld))
    ## VarianceStabilization transformation
    ## meanSdPlot(assay(vsd))
    ## dev.off()
    
    ######################################################################
    ### Pheatmap top50 DE genes
    ######################################################################
    
    #####
    # Plotting Top 50 significant DE genes with different normalization methods: 
    select <- rownames(sign.data)[1:50]
    select <- select[!is.na(select)] ## discard NA values
    
    ## Samples
    comp_name <- comp_name$category
    comp_name <- sub("\\.", "-", comp_name)
    
    numerator <- numerator$cmp1
    numerator <- sub("\\.", "-", numerator)
    
    denominator <- denominator$cmp2
    denominator <- sub("\\.", "-", denominator)
    
    print("Comparison: ")
    print(comp_name)
    
    print("numerator: ")
    print(numerator)
    
    print("denominator: ")
    print(denominator)
    
    
    print("Samples: ")
    subsheet <- df_treatment_Ind[comp_name]
    colnames(subsheet)[1] <- 'comp_name'
    #print(subsheet)
    
    listOfSampls <- c(rownames(subset(subsheet, comp_name==numerator)),
                      rownames(subset(subsheet, comp_name==denominator)))
    
    print(listOfSampls)
    
    if ( length(select) > 5 ) {
      
      ## plot rld
      try(plot1 <- pheatmap(assay(rld)[select,], main="Log Transformation Pheatmap (p.adj<0.05 and [FC]>1.2)",
                            cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                            annotation_col = df_treatment_Ind,
                            color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
                            scale="row" ## centered and scale values per row
      ))
      HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                    name_file = paste0(file_name, "_top50_DEgenes_Heatmap-LogTransformation_allSamples"), 
                                    plot_given = plot1)
      
      ## plot vsd
      #pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz.pdf")), width=15, height=12)
      try(plot2 <- pheatmap(assay(vsd)[select,], main="Variance Stabilization Pheatmap (p.adj<0.05 and [FC]>1.2)",
                            cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                            annotation_col = df_treatment_Ind,
                            color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
                            scale="row" ## centered and scale values per row7
      ))
      
      HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                    name_file = paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz_allSamples"), 
                                    plot_given = plot2)
      
      ## Only samples included in comparison
      
      ## plot rld
      try(
        
        print("select:")
        print(select)
        dataSubset <- assay(rld)[select,listOfSampls]
        print(head(dataSubset))
        
        plot3 <- pheatmap(dataSubset, main="Log Transformation Pheatmap (p.adj<0.05 and [FC]>1.2)",
                            cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                            annotation_col = df_treatment_Ind,
                            color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
                            scale="row" ## centered and scale values per row
                          )
        save_pdf(folder_path = OUTPUT_Data_sample, 
                                    name_file = paste0(file_name, "_top50_DEgenes_Heatmap-LogTransformation"), 
                                    plot_given = plot3)
      
      ## plot vsd
      #pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz.pdf")), width=15, height=12)
      
      plot4 <- pheatmap(dataSubset, main="Variance Stabilization Pheatmap (p.adj<0.05 and [FC]>1.2)",
                            cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                            annotation_col = df_treatment_Ind,
                            color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
                            scale="row" ## centered and scale values per row7
                            )
      save_pdf(folder_path = OUTPUT_Data_sample, 
                                    name_file = paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz"), 
                                    plot_given = plot4)
      )
      
    }
    
    ## Add PCA for all significant results
    data2pca <- t(sign.data[,listOfSampls])
    pca_res <- stats::prcomp(as.matrix(data2pca))
    
    pdf(file.path(OUTPUT_Data_sample,"PCA_multiple.pdf"))
    for (i in colnames(df_treatment_Ind)) {
      p<-autoplot(pca_res, 
                  data=df_treatment_Ind[listOfSampls,], colour=i) + 
        theme_classic() + ggtitle(paste0("Variable: ", i )) 
      print(p)
    }
    
    p<-autoplot(pca_res,
                data=df_treatment_Ind[listOfSampls,], colour=i) + 
      geom_text(label=listOfSampls) + 
      theme_classic() + ggtitle(paste0("Variable: ", i )) 
    print(p)
    
    dev.off()
    
    ######################################################################
    print ("Finish here for: ")
    print(file_name)
    ######################################################################
    
    data2save <- list(
      "alldata" = alldata,
      "data2pca" = data2pca,
      "rld" = rld,
      "vsd" = vsd,
      "volcan_plot" = volcan_plot,
      #"dds_object" = dds_object,
      #"res"=res,
      "res_filtered"=res_filtered,
      "sign.df"=sign.data,
      "sign.genes"=sign.data$Gene,
      "sign.count"=length(sign.data$Gene)
    )
    
    ## dump in disk RData
    save(data2save, file=file.path(OUTPUT_Data_sample, "data2return.RData"))
    
    data2return <- list(
      "res_filtered"=res_filtered,
      "sign.df"=sign.data,
      "sign.genes"=sign.data$Gene,
      "sign.count"=length(sign.data$Gene)
    )
    
  }
  
  #####
  return(data2return)  
  
  
}


#' Plot batch effect
#'
#' This functions plots original PCA and batch corrected given two variables and a putative batch variable
#' @param var1 DESeq2 object
#' @param var2 DESeq2 object
#' @param dds_object DESeq2 object
#' @param dirName Folder path to store results
#' @param batch_var Putative batch variable
#' @export
plot_batch_effect <- function(var1, var2, dds_object, dirName, batch_var) {
  ## remove batch effect?
  ## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  ## Why after VST are there still batches in the PCA plot?
  
  ## Mickael Love:
  ## DESeq() only takes as input the original counts, and this is on purpose. 
  ## This is the optimal statistical approach. To account for batch, you put 
  ## the variables at the beginning of the design, e.g. ~batch + condition 
  
  ## plot batch effect
  vsd_Test <- varianceStabilizingTransformation(dds_object, blind = FALSE)
  vsd_Test1 <- vsd_Test
  
  mat_Test <- assay(vsd_Test)
  mat_Test <- limma::removeBatchEffect(mat_Test, vsd_Test[[batch_var]])
  assay(vsd_Test) <- mat_Test
  
  p1 <- plotPCA(vsd_Test1, intgroup=var1)
  p2 <- plotPCA(vsd_Test, intgroup=var1)
  p3 <- plotPCA(vsd_Test1, intgroup=var2)
  p4 <- plotPCA(vsd_Test, intgroup=var2)
  
  pdf(paste0(dirName, "BatchEffect.pdf"))
  par(mfrow=c(1,2))
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
  ##
}

#' Adjust sample names
#'
#' Adjust samples between sample sheet and files for DESeq2
#' @param counts Expression counts. Samples as columns
#' @param target Phenotypic information. Samples as rows
#' @export
adjust_samples <- function(counts, target){
  
  ## adjust samples in data and in sample sheet
  counts <- counts[, sort(colnames(counts)) ]
  counts <- counts[, colnames(counts) %in% rownames(target) ]
  
  print ("** Samples in counts:")
  print (colnames(counts))
  print ("** Samples in target:")
  print (rownames(target))
  
  print ("** Adjusting...")
  ##
  logical_vec <- rownames(target) %in% colnames(counts)
  target <- target[logical_vec,]
  counts <- counts[, rownames(target) ]
  
  print ("** Match:")
  print (match(rownames(target), colnames(counts)))
  
  list2return <- list("counts" = counts, 
                      "target" = target)
  
  return(list2return)
  
}

#' Plot values per genes
#' 
#' Plot DESeq2 normalization values using ggplot2
#' @param gene Gene ID
#' @param tableCounts Expression counts. Samples as columns
#' @param targetsFile Phenothypic information. Sampes as rows.
#' @param condition List of phenotypic condition to retrieve in targestFile
#' @param out_folder Out folder to create plot
#' @export
plot_gene_values <- function(gene, tableCounts, targetsFile, condition, out_folder) {
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  
  gene_table <- tableCounts[tableCounts$Gene==gene,]
  print(gene_table)
  
  ## check if DE table and contains Gene, pvalue, etc
  gene_table <- dplyr::select(gene_table, -c("Gene", "baseMean", "log2FoldChange", "lfcSE","pvalue","padj"))
  print(gene_table)
  
  ## melt data
  gene_table_long <- melt(gene_table,variable.name="sample", value.name="Count")
  rownames(gene_table_long) <- gene_table_long$sample
  print(gene_table_long)
  
  ## get annotation
  print(class(condition))
  print(condition)
  annot_info <- targetsFile[rownames(gene_table_long), condition]
  print(annot)
  
  test <- merge(gene_table_long, annot_info, by.x='sample', by.y="Sample_Name")
  print(test)
  
  ## check if multiple conditions
  x.var = rlang::sym(condition[1])
  y.var = rlang::sym(condition[2])
  
  p <- ggplot(test, aes(x = !! x.var, y = Count, fill= !! y.var)) + 
    geom_boxplot() + geom_point(position=position_jitterdodge(),alpha=0.3)
  
  ## save plot
  save_pdf(folder_path = out_folder, name_file = paste0('boxplot_gene-', gene), plot_given = p)
  return(p)
}

#' Get names provided in comparison
#' 
#' When running DESeq2 you usually get names from resultsNames() such as var_comp1_vs_comp2 e.g. Sex_male_vs_female. This functions returns the name of the variables and the comparison studied.
#' @param str_given A string with the comparison. E.g. Sex_male_vs_female
#' @export
get_comparison_resultsNames <- function(str_given) {
  list_produced <- unlist(strsplit(str_given, split="_")) 
  
  str2return <- list(
    "category" = NULL,
    "cmp1" = NULL,
    "cmp2" = NULL
  )
  
  if ("vs" %in% list_produced) {
    if (list_produced[3] == "vs") {
      # category comp1 vs comp2
      str2return$category = list_produced[1]
      str2return$cmp1 = list_produced[2]
      str2return$cmp2 = list_produced[4]
    } else {
      vs_index <- as.numeric(match("vs", list_produced))
      len_given <- length(list_produced[vs_index:length(list_produced)])
      if (len_given == 2) {
        str2return$category = paste0(list_produced[1], "_", list_produced[2])
        str2return$cmp1 = list_produced[3]
        str2return$cmp2 = list_produced[5]
      } else if (len_given == 3) {
        str2return$category = list_produced[1]
        str2return$cmp1 = paste0(list_produced[2], "_", list_produced[3])
        str2return$cmp2 = paste0(list_produced[5], "_", list_produced[6])
      }
    }
  } else {
    print("Interaction term selected:")
    print(str_given)
    
    str2return$category = "interaction"
    str2return$cmp1 = str_given
    str2return$cmp2 = "reference"
    
  }
  
  return(str2return)
}


#' Relevel and rung DESEQ2 analysis
#' 
#' When running DESeq2 you sometimes require to relevel some comparisons
#' @param dds_object DESeq2 object
#' @param category Name of the variable to use within dds_object metadata
#' @param reference Name of the reference class to set.
#' @param given_dir Output dir to use
#' @param dfAnnotation Dataframe containing additional information for each sample
#' @param int_threads Number of CPUs to use [Default: 2].
#' @export
relevel_function <- function(dds_object, category, reference, given_dir, dfAnnotation, int_threads=2){
  dds_object[[category]] <- relevel(dds_object[[category]], ref=reference)
  dds_object_releveled <- DESeq(dds_object, parallel = TRUE)
  
  # check
  print("resultsNames(dds_object)")
  print(resultsNames(dds_object))
  
  print("resultsNames(dds_object_releveled)")
  print(resultsNames(dds_object_releveled))
  
  results_list <- list()
  
  for (coef_name in resultsNames(dds_object_releveled)) {
    if (coef_name=="Intercept") {} else {
      print(paste0(" + Analysis for coefficient: ", coef_name))
      listNames <- get_comparison_resultsNames(coef_name)
      
      res_dds = DESeq2_HCGB_function( 
        dds_object = dds_object_releveled, 
        coef_n = coef_name, comp_ID="relevel", comp_name = listNames[1], 
        numerator = listNames[2], denominator = listNames[3],
        OUTPUT_Data_dir = given_dir, df_treatment_Ind = dfAnnotation, 
        threads = as.numeric(int_threads))
      
      ## save to return
      #print (head(res_dds))
      results_list[[coef_name]] <- res_dds
    }
  }
  
  ## Init data to return
  data2return <- list(
    "dds_obj" = dds_object_releveled,
    "resultsNames" = resultsNames(dds_object_releveled),
    "results" = results_list
  )
  
  return(data2return)
}

#' Filter significant hits from DESEQ2 analysis
#' 
#' When running DESeq2 you usually required significant hits.
#' @param dataF Dataframe with either normalized values and DESEQ2 values or only DESEQ2.
#' @param sign_value Pvalue adjusted cutoff: Default=0.05
#' @param LFC Log Fold Change cutoff: Default: 0.26
#' @export
filter_signficant_DESEQ <- function(dataF, sign_value = 0.05, LFC=0.26) {
  
  #log2FoldChange
  #padj
  dataFilt <- subset(dataF, abs(log2FoldChange)>LFC & padj<sign_value)
  dataFilt <- dataFilt[order(dataFilt$padj),]
  return(dataFilt)
}

#' Get results from DESeq2 and normalized data
#' 
#' When running DESeq2 you usually require to get all statiscal results and normalized data
#' @param dds_obj DESeq2 object (DESeqDataSet)
#' @param coef_n Name of the coefficient obtain from DESeq2::resultsNames(dds_obj))
#' @param type By default RNA is expected and genes, either EntrezID or ENSEMBL ID is used as ID. You can specify XICRA if miRNA is provided and Gene is the combination of miRNA, variant and isomir (e.g hsa-let-7a-2-3p&iso_add3p:1&iso-23-NLJ18XQZD2). If type XICRA provided, columns is splitted into three new columns (parent, variant and UID)
#' @export
get_all_data_DESeq2 <- function(dds_obj, coef_n, type="DESeq2") {
  #res <- DESeq2::results(dds_object, name = resultsNames(dds_object)[coef_n])
  res <- DESeq2::results(dds_obj, name = coef_n)
  
  ## Merge normalized values and differential expression
  alldata <- merge(as.data.frame(counts(dds_obj, normalized=TRUE)), as.data.frame(res), by="row.names", sort=FALSE)
  names(alldata)[1] <- "Gene"
  
  if (type=="XICRA") {
    alldata <- tidyr::separate(alldata, Gene, c('parent', 'variant', 'UID'), sep = '&', remove=FALSE)  
  }
  return(alldata)
}


#' Check the effect of reducing one term from a DESEQ2 design
#' 
#' Check reducing effect using LRT method.
#' @param dds_obj.given DESeq2 object (DESeqDataSet)
#' @param formula_given Formula to substract for formula in dds_obj.given
#' @param comp.folder.given Absolute path to store results
#' @param compID.given Tag name to include for each comparison
#' @param dfAnnotation.given Dataframe with useful metadata to include
#' @export
check_reduced_LRT <- function(dds_obj.given, formula_given, 
                              comp.folder.given, compID.given, dfAnnotation.given, int_threads=2) {
  
  ## LRT: Check reduction
  DEseq.red <- DESeq2::DESeq(object = dds_obj.given, test="LRT", 
                             reduced=as.formula(formula_given))
  
  print(length(resultsNames(DEseq.red)))
  term2use <- tail(resultsNames(DEseq.red), 1)
  print(term2use)
  
  Resultsnames2use <- HCGB.IGTP.DAnalysis::get_comparison_resultsNames(term2use)
  print(Resultsnames2use)
  
  ## get for all 
  DEseq.red.res <- get_Results_DDS(dds_object = DEseq.red, 
                                   OUTPUT_Data_dir_given = comp.folder.given, 
                                   comp_ID = compID.given,
                                   dfAnnotation = dfAnnotation.given, int_threads = int_threads)
  
  return(DEseq.red.res)
}

#' Check the effect of variables in matrix design
#' 
#' When running DESeq2 you usually add multiple terms to the matrix design. Test the effect of them
#' @param sampleSheet.given Samplesheet containing metadata information
#' @param countsGiven Dataframe/matrix of counts
#' @param list.terms List of terms from samplesheet to include in design matrix
#' @param red.formula.given Design formula to use as naive and in reduction
#' @param compID.given Tag name to include for each comparison
#' @param dfAnnotation.given Dataframe with useful metadata to include
#' @param comp.folder.given Absolute path to store results
#' @param int_threads Number of threads to use

#' @export
check_terms_matrix <- function(sampleSheet.given, countsGiven, list.terms, red.formula.given, 
                               compID.given, dfAnnotation.given, comp.folder.given, int_threads=2) {
  
  resulst_list <- list()
  
  ##--------------------------
  ## naive
  ##--------------------------
  print("++++++++++++++++++++++++++++")
  print("Naive")
  print("++++++++++++++++++++++++++++")
  print("Formula: ")
  print(paste0("~", red.formula.given))
  
  naive_res <- analysis_DESeq(OUTPUT_Data_dir_given = comp.folder.given,
                              count_table = countsGiven, sample_sheet_given = sampleSheet.given, 
                              int_threads = int_threads, dfAnnotation = dfAnnotation.given, 
                              formula_given = as.formula(paste0("~", red.formula.given)), 
                              early_return = FALSE, comp_ID = paste0(compID.given, ".naive"))
  
  resulst_list[[red.formula.given]] = naive_res
  
  print("##################")
  ##--------------------------
  
  ##--------------------------
  ## add terms: addition
  ##--------------------------
  print("++++++++++++++++++++++++++++")
  print("Test addition")
  print("++++++++++++++++++++++++++++")
  
  for (term in list.terms) {
    print("Testing the effect of:")
    print(term)
    print("Formula: ")
    print(paste0("~", term, "+", red.formula.given))
    
    comp_ID.here <- paste0(compID.given, ".", term)
    
    this.term.res <- analysis_DESeq(OUTPUT_Data_dir_given = comp.folder.given,
                                    count_table = countsGiven, 
                                    sample_sheet_given = sampleSheet.given, 
                                    int_threads = int_threads, 
                                    dfAnnotation = dfAnnotation.given, 
                                    formula_given = as.formula(paste0("~", term, "+", red.formula.given)), 
                                    early_return = FALSE, 
                                    comp_ID = comp_ID.here)
    
    resulst_list[[term]] = this.term.res
    
    print("##################")
    
    print("Testing the effect of reducing:")
    print(term)
    print("Formula: ")
    print(paste0("~", term, "+", red.formula.given, " vs. ~", red.formula.given))
    
    
    ## Check reduction
    
    resulst_list[[paste0(term, ".red")]] = check_reduced_LRT(dds_obj.given = this.term.res$dds_obj, 
                                                             formula_given = as.formula(paste0("~", red.formula.given)),
                                                             comp.folder.given = comp.folder.given, 
                                                             compID.given = paste0(comp_ID.here, ".red"), 
                                                             dfAnnotation.given = dfAnnotation.given)
    
  }
  ##--------------------------
  
  ##--------------------------
  ## complex
  ##--------------------------
  ##
  print("++++++++++++++++++++++++++++")
  print("Complex:")
  print("++++++++++++++++++++++++++++")
  
  print("Formula: ")
  print(paste0("~", paste(list.terms, "+ ", collapse = ""), red.formula.given))
  
  resulst_list[["complex"]] = analysis_DESeq(OUTPUT_Data_dir_given = comp.folder.given,
                                             count_table = countsGiven, sample_sheet_given = sampleSheet.given, 
                                             int_threads = int_threads, dfAnnotation = dfAnnotation.given, 
                                             formula_given = as.formula(paste0("~", paste(list.terms, "+ ", collapse = ""), red.formula.given)), 
                                             early_return = FALSE, 
                                             comp_ID = paste0(compID.given, ".complex"))
  
  ##--------------------------
  
  print("##################")
  
  ##
  
  return(resulst_list)
}


#' DESEQ2 analysis pipeline
#' 
#' When running DESeq2 you usually add multiple terms to the matrix design. Test the effect of them
#' @param sample_sheet_given Samplesheet containing metadata information
#' @param count_table Dataframe/matrix of counts
#' @param OUTPUT_Data_dir_given Absolute path to store results
#' @param dfAnnotation Dataframe with useful metadata to include
#' @param int_threads Number of threads to use
#' @param formula_given Design formula to use
#' @param coef_n Number of the coefficient of results to test [if desired]
#' @param early_return Whether to return exploratory results early or not
#' @param comp_ID Tag name to include for each comparison
#' @param cutoff.given add an option to include cutoff when removing Zeros
#' @export
analysis_DESeq <- function(OUTPUT_Data_dir_given, count_table, sample_sheet_given, 
                           dfAnnotation, formula_given, int_threads=2,
                           coef_n=NA, early_return=FALSE, comp_ID=NULL, cutoff.given=0.9) {
  
  dir.create(OUTPUT_Data_dir_given, showWarnings = FALSE)
  
  #############
  ## Create list object for DESeq: remove 0 values
  #############
  data_DESeq <- list(
    "counts"=HCGB.IGTP.DAnalysis::discard_0_counts(countsF = count_table, cutoff = cutoff.given),
    "target"=sample_sheet_given
  )
  data_DESeq <- HCGB.IGTP.DAnalysis::adjust_samples(data_DESeq$counts, data_DESeq$target)
  
  
  ## Set parallel threads
  print (paste0("Set Multicore: ", int_threads))
  register(MulticoreParam(int_threads))
  #############
  
  #############
  ## Design ###
  #############
  ddsFullCountTable <- DESeq2::DESeqDataSetFromMatrix(
    countData = data_DESeq$counts,
    colData = data_DESeq$target, design = as.formula(formula_given) )
  
  dds_object <- DESeq2::DESeq(ddsFullCountTable, parallel = TRUE)
  
  
  ## check
  print("resultsNames(dds_object)")
  print(resultsNames(dds_object))
  
  #############
  ## exploratory dds_object
  #############
  print("Exploratory plots")
  exploratory_plots_dir <- file.path(OUTPUT_Data_dir_given, paste0(comp_ID, "_exploratory_plots"))
  dir.create(exploratory_plots_dir, showWarnings = FALSE)
  
  print(dim(sample_sheet_given))
  
  exploratory_plots_returned <- exploratory_plots(dds_object.exp = dds_object, 
                                                  OUTPUT_dir = exploratory_plots_dir, 
                                                  dfAnnotation_df = sample_sheet_given, 
                                                  list_of_cols = colnames(dfAnnotation))
  print('Out Exploratory plots')
  #############
  
  if (early_return) {
    return(list("dds_object"=dds_object, 
                "exploratory_plots" = exploratory_plots_returned,
                "resultsNames" = resultsNames(dds_object)
    ))
  }
  #############
  
  ###########
  # Get results
  #############
  results_list = get_Results_DDS(dds_object = dds_object, 
                                 OUTPUT_Data_dir_given = OUTPUT_Data_dir_given, 
                                 dfAnnotation = dfAnnotation, comp_ID = comp_ID, 
                                 int_threads = int_threads, coef_n = coef_n)
  #############
  
  #############
  ## Init data to return
  #############
  data2return <- list(
    "dds_obj" = dds_object,
    "exploratory_plots" = exploratory_plots_returned,
    "resultsNames" = resultsNames(dds_object),
    "dataDESeq" = data_DESeq,
    "formula" = formula_given,
    "results" = results_list
  )
  
  return(data2return)
}

#' Exploratory plots for DESEQ2 analysis
#' 
#' When running DESeq2 you usually add multiple terms to the matrix design. Test the effect of them
#' @param dds_object.exp DESeq2 object (DESeqDataSet)
#' @param OUTPUT_dir Absolute path to store results
#' @param dfAnnotation_df Dataframe with useful metadata to include
#' @param list_of_cols List of columns of interest to subset from metadata
#' @export
exploratory_plots <- function(dds_object.exp, OUTPUT_dir, dfAnnotation_df, list_of_cols){
  
  library(reshape2)
  library(RColorBrewer)
  
  ############################
  # Exploratory plots 
  ############################
  
  print('Inside exploratory plots')
  
  # Dispersion plot 
  plotDisp <- plotDispEsts(dds_object.exp)
  jpeg(file.path(OUTPUT_dir, "general_dispersion_plot.jpeg"))
  print(plotDisp)
  dev.off()
  
  # Top 50 genes:
  select <- order(rowMeans(counts(dds_object.exp,normalized=TRUE)),decreasing=TRUE)[1:50]
  vsd <- varianceStabilizingTransformation(dds_object.exp, blind = FALSE)
  
  top50heatmap <- pheatmap(assay(vsd)[select,], annotation = dfAnnotation_df[,list_of_cols])
  save_pdf(folder_path = OUTPUT_dir, name_file = "top50_heatmap",
           plot_given = top50heatmap)
  
  ## ape library
  # Clustering:
  # Get sample-to-sample distances
  distsRL <- dist(t(assay(vsd)))
  hc <- hclust(distsRL, method="average")
  
  # Clustering:
  ## ape library
  #pdf(file.path(OUTPUT_dir,"clustering.pdf"))
  #phylo_plot <- plot.phylo(as.phylo(hc), 
  #                         tip.color = HCGB.IGTP.DAnalysis::create_col_palette(dfAnnotation_df$condition, 
  #                                                                             levels(dfAnnotation_df$condition), 
  #                                                                             palette_given = "Set1"),
  #                         direction = "downwards", 
  #                         srt = 180, adj = 1, 
  #                         main=paste("Correlation-based clustering"),)
  #dev.off()
  
  # Get sample-to-sample distances
  distsRL <- dist(t(assay(vsd)))
  mat <- as.matrix(distsRL)
  sampleDist <- pheatmap(mat, annotation = dfAnnotation_df[,list_of_cols])
  save_pdf(folder_path = OUTPUT_dir, name_file = "heatmap_samplesDistance",
           plot_given = sampleDist)
  # 
  # ## PCA
  pdf(file.path(OUTPUT_dir,"PCA_multiple.pdf"), paper = "A4r", width = 35, height = 12)
  list_pca <- list()
  for (gr in list_of_cols) {
    plt_pca <- plotPCA(vsd, intgroup=gr) + ggtitle(gr) + 
      ggrepel::geom_text_repel(label=rownames(dfAnnotation_df)) + theme_classic()
    list_pca[[gr]] <- plt_pca
    print(plt_pca)
  }
  dev.off()
  # 
  
  #PCA_data <- plotPCA(vsd, returnData=TRUE)
  
  ### cooks distance
  df.cooks <- as.data.frame(log10(assays(dds_object.exp)[["cooks"]])) %>% melt()
  
  ## return
  plots2return <- list(
    "plotDisp" = plotDisp,
    "top50heatmap" = top50heatmap,
    #"phylo_plot" = phylo_plot,
    "sampleDist" = sampleDist,
    "PCA" = list(
      "PCA_data" = vsd,
      "PCA_list" = list_pca
    ),
    "cooks.data" = df.cooks
  )
  
  return(plots2return)
  
}


#' Generate results given a DDS object
#' 
#' When running DESeq2 you usually add multiple terms to the matrix design. Test the effect of them
#' @param dds_object DESeq2 object (DESeqDataSet)
#' @param OUTPUT_Data_dir_given Absolute path to store results
#' @param dfAnnotation Dataframe with useful metadata to include
#' @param comp_ID Tag name to include for each comparison
#' @param int_threads Number of threads to use in the analysis
#' @param coef_n Number of the coefficient of results to test [if desired]
#' @export
get_Results_DDS <- function(dds_object, OUTPUT_Data_dir_given, dfAnnotation, comp_ID,
                            int_threads=2, coef_n=NA) {
  ###########
  # Get results
  #############
  
  results_list <- list()
  if (is.numeric(coef_n)) {
    listNames <- get_comparison_resultsNames(resultsNames(dds_object)[coef_n])
    print(listNames)
    print(paste0(" + Analysis for coefficient given: ", as.character(coef_n)))
    res_dds = DESeq2_HCGB_function(
      dds_object = dds_object, coef_n = coef_n, comp_name = listNames[1], comp_ID = comp_ID,
      numerator = listNames[2], denominator = listNames[3],
      OUTPUT_Data_dir = OUTPUT_Data_dir_given, df_treatment_Ind = dfAnnotation, 
      threads = as.numeric(int_threads))
    
    ## save to return
    coef_name = as.character(resultsNames(dds_object)[coef_n])
    results_list[[coef_name]] = res_dds
    
  } else {
    
    for (coef_name in resultsNames(dds_object)) {
      if (coef_name=="Intercept") {} else {
        print(paste0(" + Analysis for coefficient: ", coef_name))
        listNames <- get_comparison_resultsNames(coef_name)
        
        print(listNames)
        
        res_dds = DESeq2_HCGB_function(
          dds_object = dds_object, coef_n = coef_name, comp_ID = comp_ID,
          comp_name = listNames[1], numerator = listNames[2], denominator = listNames[3],
          OUTPUT_Data_dir = OUTPUT_Data_dir_given, df_treatment_Ind = dfAnnotation, 
          threads = as.numeric(int_threads))
        
        ## save results
        results_list[[coef_name]] = res_dds
      }
    }
    
  }
  #############
  
  return(results_list)
}

#' Check the rank of a design matrix
#' 
#' When running DESeq2 you need a design matrix, check the rank of it first
#' @param formula2test String with formula to check
#' @param data.df Sample sheet dataframe
#' @export
check_rank_design <- function(formula2test, data.df) {
  m <- model.matrix(as.formula(formula2test), data=data.df)
  
  print("colnames(m)")
  print(colnames(m))
  
  print("Check if colSums or rowSums == 0")
  
  print("## check rows: samples")
  row.res <- apply(m, 1, function(x) all(x==0)) ## check rows: samples
  print(table(row.res))
  print("")
  print(row.res)
  print("")
  print(which(row.res))
  print("which(row.res)")
  
  print("")
  print("## check columns: categories")
  col.res <- apply(m, 2, function(x) all(x==0)) ## check columns: categories
  
  print("table(col.res)")
  print(table(col.res))
  print("")
  print(col.res)
  print("")
  print("which(col.res)")
  print(which(col.res))
}
