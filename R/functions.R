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

#' Discard low count in a dds object
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
#' >cutoff (%) of 0's. Default cutoff value is 0.9
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
#' @param threads Number of CPUs to use [Default 2].
#' @export
DESeq2_HCGB_function = function(dds_object, coef_n, name, 
                       numerator="example1", denominator="example2", 
                       OUTPUT_Data_dir, df_treatment_Ind, threads=2) {
  
  ## set name
  file_name <- paste0(name, "_", numerator, "_vs_", denominator)
  
  ## start
  print (paste0("## Starting: ", file_name))
  
  ## Set parallel threads
  print (paste0("Set Multicore: ", as.numeric(threads)))
  register(MulticoreParam(as.numeric(threads)))
  
  ## Create folder
  OUTPUT_Data_sample = file.path(OUTPUT_Data_dir, file_name)
  print (paste0("Create folder: ", OUTPUT_Data_sample))
  dir.create(OUTPUT_Data_sample)
  
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
  Pvalues_pdf <- file.path(OUTPUT_Data_sample, paste0(name,"_", numerator, "_vs_", denominator," p-value_distribution.pdf"))
  HCGB.IGTP.DAnalysis::plot_DESeq2_pvalues(Pvalues_pdf, res, res_filtered)
  
  ######################################################################
  # Write Results tables
  ######################################################################
  write.table(res_filtered, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_table.txt")), sep="\t", row.names=T, col.names=NA, quote=F)
  
  # Results table ordered by adjusted p-value
  resOrdered <- res_filtered[order(res_filtered$padj),]
  write.table(resOrdered, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_padj-ordered_table.txt")), sep="\t", row.names=T, col.names=NA, quote=F)
  
  # Save normalized values
  # Normalized values
  normValues <- counts(dds_object, normalized=T)
  write.table(normValues, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_NormValues_table.txt")), sep="\t", row.names=T, col.names=T, quote=F)
  
  ## Merge normalized values and differential expression
  alldata <- merge(as.data.frame(counts(dds_object, normalized=TRUE)), as.data.frame(res_filtered), by="row.names", sort=FALSE)
  names(alldata)[1] <- "Gene"
  ## Write results
  write.table(alldata, file=file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_NormValues-and-DE_table.txt")), sep="\t", row.names=T, col.names=NA, quote=F)
  
  ##get significant data only with counts in all samples
  ##get significant data only padj < 0.05 and log2FoldChange >1.2
  sign.data<-alldata[alldata$padj<0.05,]
  sign.data<-sign.data[!is.na(sign.data$padj),]
  sign.data<-sign.data[abs(  sign.data$log2FoldChange)>log2(1.2),]
  sign.data <- sign.data[order(sign.data$padj),]
  row.names(sign.data) <- sign.data$Gene
  
  # Results table in the same order than counting table
  write.table(sign.data, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_table_SignificantDE.txt")), sep="\t", row.names=T, col.names=NA, quote=F)
  
  ######################################################################
  ## ma plot
  ######################################################################
  jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-maplot.jpeg")), 1500, 1000, pointsize=20)
  plot_main_title <- paste0("MA Plot: ", numerator, " vs ", denominator)
  HCGB.IGTP.DAnalysis::maplot(res = alldata, main=plot_main_title)
  dev.off()
  
  ######################################################################
  ## volcano
  ######################################################################
  jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-volcano-plot.jpeg")), 1500, 1000, pointsize=20)
  volcano_main_title <- paste0("Volcano Plot: ", numerator, " vs ", denominator)
  HCGB.IGTP.DAnalysis::volcanoplot(res = alldata, main=volcano_main_title ,lfcthresh=round(log2(1.2),2), sigthresh=0.05, textcx=.8)
  dev.off()
  
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
  
  ntd <- normTransform(dds_object)
  jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-NormalTransformation-plot.jpeg")), 1500, 1000, pointsize=20)
  meanSdPlot(assay(ntd))
  dev.off()
  
  jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-RegularizedLogTransform-plot.jpeg")), 1500, 1000, pointsize=20)
  meanSdPlot(assay(rld))
  dev.off()
  
  jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-VarianceStabilizTransform-plot.jpeg")), 1500, 1000, pointsize=20)
  meanSdPlot(assay(vsd))
  dev.off()
  
  ######################################################################
  ### Pheatmap sample distribution
  ######################################################################
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
  
  pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_SampleDist.pdf")), width=15, height=12)
  pheatmap(sampleDistMatrix, 
           clustering_distance_rows=sampleDists, 
           clustering_distance_cols=sampleDists, 
           col=colors, 
           annotation_row = df_treatment_Ind)
  dev.off()
  
  ######################################################################
  ### Pheatmap top50 DE genes
  ######################################################################
  
  #####
  # Plotting Top 50 significant DE genes with different normalization methods: 
  select <- rownames(sign.data)[1:50]
  select <- select[!is.na(select)] ## discard NA values
  
  ## plot rld
  pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_top50_DEgenes_Heatmap-LogTransformation.pdf")), width=15, height=12)
  try(pheatmap(assay(rld)[select,], main="Log Transformation Pheatmap (p.adj<0.05 and [logFC]>1.2)",
	   cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
	   annotation_col = df_treatment_Ind,
	   color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
	   scale="row" ## centered and scale values per row
  ))
  dev.off()
  
  ## plot vsd
  pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz.pdf")), width=15, height=12)
  try(pheatmap(assay(vsd)[select,], main="Variance Stabilization Pheatmap (p.adj<0.05 and [logFC]>1.2)",
	   cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
	   annotation_col = df_treatment_Ind,
	   color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
	   scale="row" ## centered and scale values per row7
  ))
  dev.off()
  
  ######################################################################
  print ("Finish here for: ")
  print(file_name)
  ######################################################################
  
  #####
  return(res_filtered)
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

#' Adjust samples between sample sheet and files for DESeq2
#'
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
  
  list2return <- c("counts" = counts, 
                   "target" = target)
  
  return(list2return)
  
}
