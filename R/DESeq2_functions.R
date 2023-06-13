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

#' Plot DESeq2 p-values
#'
#' This functions plots pvalues original and after filtering
#' @param dds_object DESeq2 object
#' @param coef_n Coefficient number obtain using resultsName(dds)
#' @param comp_name Name of the comparison
#' @param comp_ID ID tag of the comparison
#' @param numerator Name of the numerator comparison
#' @param denominator Name of the denominator comparison
#' @param OUTPUT_Data_dir Folder path to store results
#' @param df_treatment_Ind Dataframe containing additional information for each sample
#' @param list_of_cols Set of columns with important information in df_treatment_ind 
#' @param threads Number of CPUs to use [Default: 2].
#' @param sign_value.given Adjusted pvlaue cutoff. Default=0.05, 
#' @param LFC.given Log Fold change cutoff. Default=log2(1.2), 
#' @param forceResults Boolean to force re-run analysis if already generated in the folder provided
#' @param gene.annot.df Dataframe containing gene annotation (Default: NULL)
#' @export

DESeq2_HCGB_function = function(dds_object, coef_n, comp_name, comp_ID="comp1",
                                numerator="example1", denominator="example2", 
                                OUTPUT_Data_dir, df_treatment_Ind, list_of_cols, threads=2, 
                                sign_value.given = 0.05, LFC.given = log2(1.2), gene.annot.df=NULL,
                                forceResults=FALSE, shrinkage="apeglm") {
  
  #--------------------------
  # Packages
  #--------------------------
  library(DESeq2)
  library(vsn)
  library(EnhancedVolcano)
  library(BiocParallel)
  library(ggfortify)
  library(ggrepel)
  library(pheatmap)
  library(reshape2)
  library(RColorBrewer)
  #--------------------------
  
  #--------------------------
  ## Prepare
  #--------------------------
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
      return(file.path(OUTPUT_Data_sample, "data2return.RData"))
    }
  }
  #--------------------------
  
  #--------------------------
  ## Generate results according to comparison
  #--------------------------
  
  print("+ Using shrinkage estimation provided:")
  print(shrinkage)
  
  ## Compare pvalues distribution
  res <- DESeq2::lfcShrink(dds_object, coef=coef_n, type = shrinkage, parallel = TRUE)
  
  ## filter low counts
  dds_object = HCGB.IGTP.DAnalysis::discard_lowCounts(dds_object = dds_object)
  res_filtered <- DESeq2::lfcShrink(dds_object, coef=coef_n, type=shrinkage, parallel = TRUE)
  
  #--------------------------
  
  #--------------------------
  ### Plot p-values
  #--------------------------
  ## original/filtered
  Pvalues_pdf <- file.path(OUTPUT_Data_sample, paste0(comp_name,"_", numerator, "_vs_", denominator," p-value_distribution.pdf"))
  HCGB.IGTP.DAnalysis::plot_DESeq2_pvalues(Pvalues_pdf, res, res_filtered)
  #--------------------------
  
  #--------------------------
  # Write Results tables
  #--------------------------
  
  ## Merge normalized values and differential expression
  alldata <- merge(as.data.frame(counts(dds_object, normalized=TRUE)), as.data.frame(res_filtered), by="row.names")
  rownames(alldata) <- alldata$Row.names
  alldata$Row.names <- NULL
  #names(alldata)[1] <- "Gene"
  
  if (!is.null(gene.annot.df)) {
    alldata <- merge(gene.annot.df, alldata, by='row.names')
    rownames(alldata) <- alldata$Row.names
    alldata$Row.names <- NULL
  }
  
  Gene <- rownames(alldata)
  alldata <- cbind(Gene, alldata)
  print(head(alldata))
  
  ## Write results
  write.table(alldata, file=file.path(OUTPUT_Data_sample, 
                                      paste0(file_name, "-ResultsCounting_NormValues-and-DE_table.txt")), 
              sep="\t", quote=T, row.names = F)
  
  ##get significant data only with counts in all samples
  sign.data <- filter_signficant_DESEQ(alldata, sign_value = sign_value.given, LFC = LFC.given)
  row.names(sign.data) <- sign.data$Gene
  
  # Results table in the same order than counting table
  write.table(sign.data, 
              file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_table_SignificantDE.txt")), 
              sep="\t", quote=T, row.names = F)
  #--------------------------
  
  #--------------------------
  ## Samples
  #--------------------------
  comp_name <- comp_name$category
  #comp_name <- sub("\\.", "-", comp_name)
  
  numerator <- numerator$cmp1
  #numerator <- sub("\\.", "-", numerator)
  
  denominator <- denominator$cmp2
  #denominator <- sub("\\.", "-", denominator)
  
  print("Hi there!")
  alldata2 <- alldata
  
  listOfSampls <- c()
  
  if (denominator=="reference") {
    
    listOfSampls = NULL
    Samplslist <- list()
    
  } else {
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
    
    print("List of samples in these comparison")
    print(listOfSampls)
    
    Samplslist <- list(
      numerator = rownames(subset(subsheet, comp_name==numerator)),
      denominator = rownames(subset(subsheet, comp_name==denominator))
    )
    
    ## Get data for this samples
    ## 
    
    ## add basemean for each category
    alldata2['baseMean_num'] <- rowMeans(alldata[,Samplslist$numerator])
    alldata2['baseMean_den'] <- rowMeans(alldata[,Samplslist$denominator])
    
    ## number of 0s
    alldata2['0counts_num'] <- apply(alldata[,Samplslist$numerator], 
                                     1, function(x) sum(x == 0))
    alldata2['0counts_den'] <- apply(alldata[,Samplslist$denominator], 
                                     1, function(x) sum(x == 0))
    
    ## percentage of counts
    alldata2['0counts.perc_num'] <- round(apply(alldata2[,Samplslist$numerator], 
                                                1, function(x) sum(x == 0))/length(Samplslist$numerator)*100, digits = 2)
    
    alldata2['0counts.perc_den'] <- round(apply(alldata2[,Samplslist$denominator], 
                                                1, function(x) sum(x == 0))/length(Samplslist$denominator)*100, digits = 2)
    
    if (!is.null(gene.annot.df)) {
      
      list_cols <- c("Gene", "ensembl_gene_id", "hgnc_symbol", "description", "gene_biotype", "baseMean", "baseMean_num", 
                     "baseMean_den", "0counts_num", "0counts.perc_num", 
                     "0counts_den", "0counts.perc_den", 
                     "log2FoldChange", "lfcSE", "pvalue", "padj")
      
    } else {
      list_cols <- c("Gene", "baseMean", "baseMean_num", 
                     "baseMean_den", "0counts_num", "0counts.perc_num", 
                     "0counts_den", "0counts.perc_den", 
                     "log2FoldChange", "lfcSE", "pvalue", "padj")
      
    }
    
    alldata2 <- alldata2[,c(list_cols, Samplslist$numerator, Samplslist$denominator)]
    
    
    library(openxlsx)
    
    DE.filename <- file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting.xlsx"))
    sheet_name <- paste0(numerator, "_vs_", denominator)
    title_name <- paste0("Comparison for: ", comp_name, ": ", numerator, " vs. ", denominator)
    
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, title_name, startRow = 2, startCol=2,rowNames = FALSE, keepNA=TRUE,na.string="NA")
    openxlsx::writeData(wb, sheet_name, alldata2, startRow = 4, startCol=2, rowNames = FALSE, keepNA=TRUE,na.string="NA")
    
    sheet_name2 <- "all.data"
    openxlsx::addWorksheet(wb, sheet_name2)
    openxlsx::writeData(wb, sheet_name2, paste0(title_name, ": all other samples included"), 
                        startRow = 2, startCol=2,rowNames = FALSE, keepNA=TRUE,na.string="NA")
    openxlsx::writeData(wb, sheet_name2, alldata, startRow = 4, startCol=2, rowNames = FALSE, keepNA=TRUE,na.string="NA")
    
    ## add data with all counts
    ##
    openxlsx::saveWorkbook(wb, DE.filename, overwrite = TRUE)
    
  }
  
  
  #--------------------------
  ## check if it is worth to continue, avoid error if missing sign.data
  #--------------------------
  print(length(rownames(sign.data)))
  
  if (length(rownames(sign.data)) < 3) {
    
    data2save<- list(
      "alldata2" = alldata2,
      "alldata" = alldata,
      "Samplslist" = Samplslist,
      #"dds_object" = dds_object,
      #"res"=res,
      "res_filtered" = res_filtered,
      "sign.df"=sign.data,
      "sign.genes"=sign.data$Gene,
      "sign.count"=length(sign.data$Gene),
      "shrinkage.LFC"=shrinkage
    )
    
    ## dump in disk RData
    save(data2save, file=file.path(OUTPUT_Data_sample, "data2return.RData"))
    
    data2return<- list(
      "alldata2" = alldata2,
      "Samplslist" = Samplslist,
      "alldata" = alldata,
      "res_filtered" = res_filtered,
      "sign.df"=sign.data,
      "sign.genes"=sign.data$Gene,
      "sign.count"=length(sign.data$Gene),
      "shrinkage.LFC"=shrinkage
    )
    
    return(data2return)
  } 
  #--------------------------  
  
  #--------------------------
  ## ma plot
  #--------------------------
  pdf(file.path(OUTPUT_Data_sample, "DiffExpression-maplot.jpeg"))
  plotMA(res_filtered)
  dev.off()
  #--------------------------
  
  #--------------------------
  ## volcano
  #--------------------------
  #jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-volcano-plot.jpeg")), 1500, 1000, pointsize=20)
  volcano_main_title <- paste0("Volcano Plot: ", numerator, " vs ", denominator)
  volcan_plot <- EnhancedVolcano::EnhancedVolcano(alldata, x="log2FoldChange", y="padj", lab="",
                                                  pCutoff=sign_value.given, FCcutoff=LFC.given, pointSize=3, labSize=6) + 
    ggplot2::scale_x_continuous() + ggplot2::labs(title = volcano_main_title)
  
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                name_file = paste0(file_name, "_DiffExpression-volcano-plot"), plot_given = volcan_plot)
  #--------------------------
  
  #--------------------------
  ## Transform normal
  #--------------------------
  
  ## rlog transformation
  #rld <- DESeq2::rlogTransformation(dds_object)
  
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
  #--------------------------
  
  #--------------------------
  ### Pheatmap top50 DE genes
  #--------------------------
  
  print("Plotting Top 50 significant DE genes with different normalization methods")
  # Plotting Top 50 significant DE genes with different normalization methods:
  select <- rownames(sign.data)[1:50]
  select <- select[!is.na(select)] ## discard NA values
  
  
  data2pheatmap <- discard_0_counts(countsF = assay(vsd), cutoff = 0.75)
  
  if ( length(select) > 5 ) {
    
    ## plot rld
    #plot1 <- pheatmap(assay(rld)[select,],
    #                  main=paste0("Log Transformation Pheatmap (p.adj<", sign_value.given, " and [LFC]>", LFC.given),
    #                  cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
    #                  annotation_col = df_treatment_Ind[,list_of_cols],
    #                  color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)),
    #                  scale="row" ## centered and scale values per row
    #)
    #HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample,
    #                              name_file = paste0(file_name, "_top50_DEgenes_Heatmap-LogTransformation_allSamples"),
    #                              plot_given = plot1)
    
    ## plot vsd
    plot2 <- pheatmap(data2pheatmap[select,],
                      main=paste0("Variance Stabilization Pheatmap (p.adj<", sign_value.given, " and [LFC]>", LFC.given),
                      cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                      annotation_col = df_treatment_Ind[,list_of_cols],
                      color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)),
                      scale="row" ## centered and scale values per row7
    )
    
    HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample,
                                  name_file = paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz_allSamples"),
                                  plot_given = plot2)
    
    
    if (!is.null(listOfSampls)) {
      
      ## Only samples included in comparison
      dataSubset <- try(data2pheatmap[select,listOfSampls], silent = TRUE)
      
      if (exists("dataSubset")) {
        
        print("select:")
        print(select)
        print(head(dataSubset))
        
        ## plot rld
        #plot3 <- pheatmap(dataSubset,
        #                  main=paste0("Log Transformation Pheatmap (p.adj<", sign_value.given, " and [LFC]>", LFC.given),
        #                  cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
        #                  annotation_col = df_treatment_Ind[,list_of_cols],
        #                  color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)),
        #                  scale="row" ## centered and scale values per row
        #)
        #save_pdf(folder_path = OUTPUT_Data_sample,
        #         name_file = paste0(file_name, "_top50_DEgenes_Heatmap-LogTransformation"),
        #         plot_given = plot3)
        
        ## plot vsd
        plot4 <- pheatmap(dataSubset,
                          main=paste0("Variance Stabilization Pheatmap (p.adj<", sign_value.given, " and [LFC]>", LFC.given),
                          cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                          annotation_col = df_treatment_Ind[,list_of_cols],
                          color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)),
                          scale="row" ## centered and scale values per row7
        )
        save_pdf(folder_path = OUTPUT_Data_sample,
                 name_file = paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz"),
                 plot_given = plot4)
      }
      
    }
    
  }
  #--------------------------
  
  #--------------------------
  ## Add PCA for all significant results
  #--------------------------
  ## create PCA
  ## Add PCA for all significant results
  norm.counts <- as.data.frame(counts(dds_object, normalized=TRUE))
  data2pca <- t(norm.counts[sign.data$Gene,])
  pca_res <- stats::prcomp(as.matrix(data2pca))
  
  pdf(file.path(OUTPUT_Data_sample,"PCA_multiple.pdf"))
  for (i in colnames(df_treatment_Ind[,list_of_cols])) {
    p<-autoplot(pca_res, 
                data=df_treatment_Ind, 
                colour=i) + 
      theme_classic() + 
      ggtitle(paste0("Variable: ", i )) 
    print(p)
  }
  
  p<-autoplot(pca_res,
              data=df_treatment_Ind, 
              colour=i) + 
    geom_text(label=rownames(df_treatment_Ind)) + 
    theme_classic() + ggtitle(paste0("Variable: ", i )) 
  print(p)
  
  dev.off()
  #--------------------------
  
  #--------------------------
  ## Add boxplot for each DE gene
  #--------------------------
  boxplot_DE <- file.path(OUTPUT_Data_sample, "boxplot_DE")
  dir.create(boxplot_DE)
  
  df_treatment_Ind
  
  DE_plots.df <- data.frame(row.names = rownames(df_treatment_Ind), 
                            df_treatment_Ind[,list_of_cols],
                            t(sign.data[,rownames(df_treatment_Ind)]))
  print(DE_plots.df)
  
  ## print only top50
  for (gene_given in head(rownames(sign.data), n=50)) {
    
    ## 
    print(gene_given)
    g <- gsub("-", "\\.", gene_given)
    g <- gsub("&", "\\.", g)
    g <- gsub(":", "\\.", g)
    g <- gsub("\\+", "\\.", g)
    
    ##
    print(g)
    if (!is.null(gene.annot.df)) {
      gene_annot.df <- gene.annot.df.clean[gene_given,]
      gene_name = paste0(gene_given, "_", gene_annot.df$hgnc_symbol)
      print(gene_name)
    } else {
      gene_name = g
    }
      
    pdf(file.path(boxplot_DE, paste0(gene_name, ".pdf")), paper = "A4r", width = 35, height = 12)
    for (i in colnames(df_treatment_Ind[,list_of_cols])) {
      g <- gsub("-", "\\.", g)
      
      print(g)
      
      if (is.numeric(df_treatment_Ind[,i])) {
        p2 <- ggscatter_plotRegression(data_all_given = DE_plots.df, x.given = g, y.given = i, title_string = i) 
      } else {
        p2 <- ggboxplot_scatter(data_all_given = DE_plots.df, colName = i, y.coord = g)   
      }
      
      print(p2)
    }
    dev.off()
    
  }
  #--------------------------
  
  ######################################################################
  print ("Finish here for: ")
  print(file_name)
  ######################################################################
  
  ######################################################################
  ## Save
  ######################################################################
  data2save <- list(
    "alldata2" = alldata2,
    "Samplslist" = Samplslist,
    "alldata" = alldata,
    "data2pca" = data2pca,
    #"rld" = rld, ## It takes too much time
    "vsd" = vsd,
    "volcan_plot" = volcan_plot,
    #"dds_object" = dds_object,
    #"res"=res,
    "res_filtered"=res_filtered,
    "sign.df"=sign.data,
    "sign.genes"=sign.data$Gene,
    "sign.count"=length(sign.data$Gene),
    "shrinkage.LFC"=shrinkage,
    "DE_plots.df" = DE_plots.df
  )
  
  ## dump in disk RData
  save(data2save, file=file.path(OUTPUT_Data_sample, "data2return.RData"))
  
  data2return <- list(
    "alldata2" = alldata2,
    "Samplslist" = Samplslist,
    "alldata" = alldata,
    "res_filtered"=res_filtered,
    "sign.df"=sign.data,
    "sign.genes"=sign.data$Gene,
    "sign.count"=length(sign.data$Gene),
    "shrinkage.LFC"=shrinkage,
    "DE_plots.df"=DE_plots.df
  )
  ######################################################################
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
  print (rownames(target) == colnames(counts))
  
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
plot_gene_values <- function(gene, tableCounts, targetsFile, condition, out_folder, print2file=TRUE) {
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
  
  if (print2file) {
    ## save plot
    save_pdf(folder_path = out_folder, name_file = paste0('boxplot_gene-', gene), plot_given = p)
  }
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
#' @param sign_value.given Adjusted pvlaue cutoff. Default=0.05, 
#' @param LFC.given Log Fold change cutoff. Default=log2(1.2), 
#' @param forceResults Boolean to force re-run analysis if already generated in the folder provided
#' @param localFit Use a fitType=local for mean dispersion fit in DESeq2
#' @param shrinkage.given LFC shrinkage estimator provided. Available: apeglm, ashr or normal
#' @export
relevel_function <- function(dds_object, category, reference, 
                             given_dir, dfAnnotation, list_of_cols, gene.annot.df.given,
                             int_threads=2, sign_value.given = 0.05, LFC.given = log2(1.2), 
                             comp_ID.given="comp1", forceResults=FALSE, localFit=FALSE, shrinkage.given='apeglm'){
  ## relevel
  dds_object[[category]] <- relevel(dds_object[[category]], ref=reference)
  
  ## re-run dispersion fit
  if (localFit) {
    dds_object_releveled <- DESeq(dds_object, parallel = TRUE, fitType = "local")  
  } else {
    dds_object_releveled <- DESeq(dds_object, parallel = TRUE)
  }
  
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
        coef_n = coef_name, comp_ID=paste0("relevel_", comp_ID.given), 
        comp_name = listNames[1], 
        numerator = listNames[2], denominator = listNames[3],
        OUTPUT_Data_dir = given_dir, df_treatment_Ind = dfAnnotation, list_of_cols = list_of_cols,
        sign_value.given = sign_value.given, LFC.given = LFC.given,
        threads = as.numeric(int_threads), gene.annot.df = gene.annot.df.given,
        forceResults = forceResults, shrinkage = shrinkage.given)
      
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
#' @param forceResults Boolean to force re-run analysis if already generated in the folder provided
#' @param shrinkage.given LFC shrinkage estimator provided. Available: apeglm, ashr or normal
#' @export
check_reduced_LRT <- function(dds_obj.given, formula_given, 
                              comp.folder.given, compID.given, 
                              dfAnnotation.given,  list_of_cols, int_threads=2, 
                              sign_value.given=0.05, LFC.given = log2(1.2),
                              gene.annot=NULL, forceResults=FALSE, shrinkage.given="apeglm") {
  
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
                                   dfAnnotation = dfAnnotation.given, list_of_cols = list_of_cols,
                                   int_threads = int_threads, 
                                   forceResults=forceResults,   
                                   sign_value.given = sign_value.given, 
                                   LFC.given = LFC.given, 
                                   gene.annot=gene.annot, shrinkage.given=shrinkage.given)
  
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
#' @param comp.folder.given Absolute path to store results
#' @param int_threads Number of threads to use
#' @param shrinkage.given LFC shrinkage estimator provided. Available: apeglm, ashr or normal
#' @export
check_terms_matrix <- function(sampleSheet.given, countsGiven, list.terms, red.formula.given, list_of_cols,
                               compID.given, comp.folder.given, 
                               int_threads=2, gene.annot.df=NULL, shrinkage.given="apeglm") {
  
  resulst_list <- list()
  
  ##--------------------------
  ## naive
  ##--------------------------
  print("++++++++++++++++++++++++++++")
  print("Naive")
  print("++++++++++++++++++++++++++++")
  print("Formula: ")
  print(paste0("~", red.formula.given))
  
  naive_res <- HCGB.IGTP.DAnalysis::analysis_DESeq(OUTPUT_Data_dir_given = comp.folder.given,
                                                   count_table = countsGiven, 
                                                   sample_sheet_given = sampleSheet.given, 
                                                   int_threads = int_threads,  
                                                   gene.annot=gene.annot.df,
                                                   list_of_cols=list_of_cols,
                                                   formula_given = as.formula(paste0("~", red.formula.given)), 
                                                   early_return = FALSE, comp_ID = paste0(compID.given, ".naive"), shrinkage.given=shrinkage.given)
  
  ## Save only the path
  save(naive_res, file = file.path(comp.folder.given, "naive.RData"))
  naive_res <- NULL
  resulst_list[['naive']] = file.path(comp.folder.given, "naive.RData")
  
  print("##################")
  ##--------------------------
  
  ##--------------------------
  ## add terms: addition
  ##--------------------------
  print("++++++++++++++++++++++++++++")
  print("Test addition and interaction")
  print("++++++++++++++++++++++++++++")
  
  for (term in list.terms) {
    
    print("##################")
    print("Testing the effect of:")
    print(term)
    
    comp_ID.here <- paste0(compID.given, "_", term)
    
    print("Testing the effect of addittion:")
    print("Formula: ")
    print(paste0("~", term, "+", red.formula.given))
    
    this.term.res.add <- try(HCGB.IGTP.DAnalysis::analysis_DESeq(OUTPUT_Data_dir_given = comp.folder.given,
                                                                 count_table = countsGiven, 
                                                                 sample_sheet_given = sampleSheet.given, 
                                                                 int_threads = int_threads, 
                                                                 formula_given = as.formula(paste0("~", term, "+", red.formula.given)), 
                                                                 early_return = FALSE, 
                                                                 gene.annot=gene.annot.df,
                                                                 list_of_cols=list_of_cols,
                                                                 comp_ID = paste0(comp_ID.here, ".add"), shrinkage.given=shrinkage.given))
    ## Save only the path
    save(this.term.res.add, file = file.path(comp.folder.given, paste0(term, ".add.RData")))
    resulst_list[[paste0(term, ".add")]] = file.path(comp.folder.given, paste0(term, ".add.RData"))
    
    print("##################")
    print("Testing the effect of interaction:")
    print("Formula: ")
    print(paste0("~", term, ":", red.formula.given))
    
    this.term.res.int <- try(HCGB.IGTP.DAnalysis::analysis_DESeq(OUTPUT_Data_dir_given = comp.folder.given,
                                                                 count_table = countsGiven, 
                                                                 sample_sheet_given = sampleSheet.given, 
                                                                 int_threads = int_threads, 
                                                                 formula_given = as.formula(paste0("~", term, ":", red.formula.given)), 
                                                                 early_return = FALSE, 
                                                                 gene.annot=gene.annot.df,
                                                                 list_of_cols=list_of_cols,
                                                                 comp_ID = paste0(comp_ID.here, ".int"), shrinkage.given=shrinkage.given))
    ## Save only the path
    save(this.term.res.int, file = file.path(comp.folder.given, paste0(term, ".int.RData")))
    resulst_list[[paste0(term, ".int")]] = file.path(comp.folder.given, paste0(term, ".int.RData"))
    
    print("##################")
    print("Testing the effect of reducing:")
    print(term)
    print("Formula: ")
    print(paste0("~", term, ":", red.formula.given, " vs. ~", red.formula.given))
    
    ## Check reduction
    red.res = try(HCGB.IGTP.DAnalysis::check_reduced_LRT(dds_obj.given = this.term.res.int$dds_obj, 
                                                         formula_given = as.formula(paste0("~", red.formula.given)),
                                                         comp.folder.given = comp.folder.given, 
                                                         gene.annot=gene.annot.df, 
                                                         list_of_cols=list_of_cols,
                                                         compID.given = paste0(comp_ID.here, ".red"), 
                                                         dfAnnotation.given = sampleSheet.given, shrinkage.given=shrinkage.given))
    ## save the path
    save(red.res, file = file.path(comp.folder.given, paste0(term, ".red.RData")))
    resulst_list[[paste0(term, ".red")]] = file.path(comp.folder.given, paste0(term, ".red.RData"))
    
    red.res <- NULL
    this.term.res.int <- NULL
    this.term.res.add <- NULL
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
  
  complex_res = analysis_DESeq(OUTPUT_Data_dir_given = comp.folder.given,
                               count_table = countsGiven, sample_sheet_given = sampleSheet.given, 
                               int_threads = int_threads, gene.annot=gene.annot.df,
                               list_of_cols=list_of_cols,
                               formula_given = as.formula(paste0("~", 
                                                                 paste(list.terms, "+ ", 
                                                                       collapse = ""), red.formula.given)), 
                               early_return = FALSE, 
                               comp_ID = paste0(compID.given, ".complex"), shrinkage.given=shrinkage.given)
  
  ## save the path
  save(complex_res, file = file.path(comp.folder.given, "complex.RData"))
  resulst_list[["complex"]] = file.path(comp.folder.given, "complex.RData")
  
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
#' @param int_threads Number of threads to use
#' @param formula_given Design formula to use
#' @param coef_n Number of the coefficient of results to test [if desired]
#' @param early_return Whether to return exploratory results early or not
#' @param comp_ID Tag name to include for each comparison
#' @param cutoff.given add an option to include cutoff when removing Zeros
#' @param sign_value.given Adjusted pvalue cutoff. Default=0.05, 
#' @param LFC.given Log Fold change cutoff. Default=log2(1.2), 
#' @param localFit Use a fitType=local for mean dispersion fit in DESeq2
#' @param forceResults Boolean to force re-run analysis if already generated in the folder provided
#' @param gene.annot Dataframe containing gene annotation (Default: NULL)
#' @param shrinkage.given LFC shrinkage estimator provided. Available: apeglm, ashr or normal
#' @export
analysis_DESeq <- function(OUTPUT_Data_dir_given, count_table, sample_sheet_given, 
                           list_of_cols, formula_given, int_threads=2,
                           sign_value.given = 0.05, LFC.given = log2(1.2),
                           coef_n=NA, early_return=FALSE, comp_ID=NULL, cutoff.given=0.9, 
                           localFit=FALSE, forceResults=FALSE, gene.annot=NULL, shrinkage.given="apeglm") {
  
  dir.create(OUTPUT_Data_dir_given, showWarnings = FALSE)
  
  #############
  ## Create list object for DESeq: remove 0 values
  #############
  data_DESeq <- list(
    "counts"=HCGB.IGTP.DAnalysis::discard_0_counts(countsF = count_table, cutoff = cutoff.given),
    "target"=sample_sheet_given
  )
  data_DESeq <- adjust_samples(data_DESeq$counts, data_DESeq$target)
  
  
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
  
  if (localFit) {
    print("fitType = 'local'")
    dds_object <- DESeq2::DESeq(ddsFullCountTable, parallel = TRUE, fitType = "local")
  } else {
    print("fitType = 'parametric'")
    dds_object <- DESeq2::DESeq(ddsFullCountTable, parallel = TRUE)  
  }
  
  ## check
  print("resultsNames(dds_object)")
  print(resultsNames(dds_object))
  
  #############
  ## exploratory dds_object
  #############
  print("Exploratory plots")
  exploratory_plots_dir <- file.path(OUTPUT_Data_dir_given, paste0(comp_ID, "_exploratory_plots"))
  dir.create(exploratory_plots_dir, showWarnings = FALSE)
  
  exploratory_plots_returned <- exploratory_plots(dds_object.exp = dds_object, 
                                                  OUTPUT_dir = exploratory_plots_dir, 
                                                 dfAnnotation_df = sample_sheet_given, 
                                                 list_of_cols = list_of_cols)
  print('Out Exploratory plots')
  #############
  
  if (early_return) {
    return(list("dds_object"=dds_object, 
                "exploratory_plots" = exploratory_plots_returned,
                "resultsNames" = resultsNames(dds_object)
    ))
  }
  #############
  
  ############ 
  # Norm.counts
  ############
  # Save normalized values
  # Normalized values
  dds_object1 = HCGB.IGTP.DAnalysis::discard_lowCounts(dds_object = dds_object)
  
  normValues <- counts(dds_object1, normalized=T)
  print(head(normValues))
  
  if (!is.null(gene.annot)) {
    normValues <- merge(gene.annot, normValues, by='row.names')
    print(head(normValues))
    rownames(normValues) <- normValues$Row.names
    normValues$Row.names <- NULL
    print("## Add annotation")
    print(head(normValues))
  }
  
  Gene <- rownames(normValues)
  normValues <- cbind(Gene, normValues)
  write.table(normValues, file.path(OUTPUT_Data_dir_given, "NormValues_table.txt"), sep="\t", row.names=F, col.names=T, quote=T)
  
  ###########
  # Get results
  #############
  results_list = get_Results_DDS(dds_object = dds_object, 
                                 OUTPUT_Data_dir_given = OUTPUT_Data_dir_given, 
                                 dfAnnotation = sample_sheet_given, list_of_cols, comp_ID = comp_ID, 
                                 sign_value.given = sign_value.given, LFC.given = LFC.given, gene.annot=gene.annot,
                                 int_threads = int_threads, coef_n = coef_n, forceResults=forceResults, shrinkage.given=shrinkage.given)
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
    "results" = results_list,
    "normValues" = normValues
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
#' @param forceResults Boolean to force re-run analysis if already generated in the folder provided
#' @export
exploratory_plots <- function(dds_object.exp, OUTPUT_dir, dfAnnotation_df, list_of_cols, forceResults=FALSE){
  
  library(reshape2)
  library(RColorBrewer)
  
  if (forceResults) {
    
  } else {
    ## check if previously done
    if (file.exists(file.path(OUTPUT_dir, "exploratory.RData"))) {
      print("Data already available in:")
      print(OUTPUT_dir)
      return(file.path(OUTPUT_dir, "exploratory.RData"))
    }
  }
  
  ############################
  # Exploratory plots 
  ############################
  
  print('Inside exploratory plots')
  
  ## Create plot for different dispersions and return value of best fit
  
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
  for (gr in colnames(dfAnnotation_df)) {
    print(paste0("Printing PCA for ", gr))
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
  
  save(plots2return, file = file.path(OUTPUT_dir, "exploratory.RData"))
  
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
#' @param sign_value.given Adjusted pvlaue cutoff. Default=0.05, 
#' @param LFC.given Log Fold change cutoff. Default=log2(1.2), 
#' @param forceResults Boolean to force re-run analysis if already generated in the folder provided
#' @param gene.annot Dataframe containing gene annotation (Default: NULL)
#' @param shrinkage.given LFC shrinkage estimator provided. Available: apeglm, ashr or normal
#' @export
get_Results_DDS <- function(dds_object, OUTPUT_Data_dir_given, dfAnnotation, list_of_cols, comp_ID,
                            sign_value.given = 0.05, LFC.given = log2(1.2),
                            int_threads=2, coef_n=NA, forceResults=FALSE, gene.annot=NULL, shrinkage.given='apeglm') {
  ###########
  # Get results
  #############
  
  results_list <- list()
  if (is.numeric(coef_n)) {
    listNames <- get_comparison_resultsNames(resultsNames(dds_object)[coef_n])
    print(listNames)
    print(paste0(" + Analysis for coefficient given: ", as.character(coef_n)))
    res_dds = try(DESeq2_HCGB_function(
      dds_object = dds_object, coef_n = coef_n, comp_name = listNames[1], comp_ID = comp_ID,
      numerator = listNames[2], denominator = listNames[3],
      OUTPUT_Data_dir = OUTPUT_Data_dir_given, df_treatment_Ind = dfAnnotation, 
      list_of_cols = list_of_cols,
      sign_value.given = sign_value.given, LFC.given = LFC.given,
      threads = as.numeric(int_threads), forceResults=forceResults, gene.annot.df = gene.annot, shrinkage=shrinkage.given))
    
    ## save to return
    coef_name = as.character(resultsNames(dds_object)[coef_n])
    results_list[[coef_name]] = res_dds
    
  } else {
    
    for (coef_name in resultsNames(dds_object)) {
      if (coef_name=="Intercept") {} else {
        print(paste0(" + Analysis for coefficient: ", coef_name))
        listNames <- get_comparison_resultsNames(coef_name)
        
        print(listNames)
        
        res_dds = try(DESeq2_HCGB_function(
          dds_object = dds_object, coef_n = coef_name, comp_ID = comp_ID,
          comp_name = listNames[1], numerator = listNames[2], denominator = listNames[3],
          OUTPUT_Data_dir = OUTPUT_Data_dir_given, df_treatment_Ind = dfAnnotation, list_of_cols = list_of_cols,
          sign_value.given = sign_value.given, LFC.given = LFC.given,
          threads = as.numeric(int_threads), forceResults=forceResults, gene.annot=gene.annot, shrinkage=shrinkage.given))
        
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


