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
  
  
  library(DESeq2)
  library(vsn)
  library(EnhancedVolcano)
  
  
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
  #jpeg(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-volcano-plot.jpeg")), 1500, 1000, pointsize=20)
  volcano_main_title <- paste0("Volcano Plot: ", numerator, " vs ", denominator)
  volcan_plot <- EnhancedVolcano::EnhancedVolcano(alldata, x="log2FoldChange", y="padj", lab="",
                                   pCutoff=0.05, FCcutoff=1.2, pointSize=3, labSize=6) + ggplot2::scale_x_continuous()
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
  
  ntd <- normTransform(dds_object)
  pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_DiffExpression-Transformation.pdf")))
  ## Normal transformation  
  meanSdPlot(assay(ntd))
  ## RLE transformation
  meanSdPlot(assay(rld))
  ## VarianceStabilization transformation
  meanSdPlot(assay(vsd))
  dev.off()
  
  ######################################################################
  ### Pheatmap top50 DE genes
  ######################################################################
  
  #####
  # Plotting Top 50 significant DE genes with different normalization methods: 
  select <- rownames(sign.data)[1:50]
  select <- select[!is.na(select)] ## discard NA values

  if ( length(select) > 5 ) {
  
  ## plot rld
  try(plot1 <- pheatmap(assay(rld)[select,], main="Log Transformation Pheatmap (p.adj<0.05 and [logFC]>1.2)",
	   cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
	   annotation_col = df_treatment_Ind,
	   color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
	   scale="row" ## centered and scale values per row
  ))
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                name_file = paste0(file_name, "_top50_DEgenes_Heatmap-LogTransformation"), plot_given = plot1)
  
  ## plot vsd
  #pdf(file.path(OUTPUT_Data_sample, paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz.pdf")), width=15, height=12)
  try(plot2 <- pheatmap(assay(vsd)[select,], main="Variance Stabilization Pheatmap (p.adj<0.05 and [logFC]>1.2)",
	   cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
	   annotation_col = df_treatment_Ind,
	   color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
	   scale="row" ## centered and scale values per row7
  ))
  
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                name_file = paste0(file_name, "_top50_DEgenes_Heatmap-VarianceStabiliz"), plot_given = plot2)
  
  }

  ######################################################################
  print ("Finish here for: ")
  print(file_name)
  ######################################################################
  
  #####
  return(alldata)
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
#' @param condition List oif phenotypic condition to retrieve in targestFile
#' @param out_folder Outfolder to create plot
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
  
  pdf(file.path(out_folder, paste0('boxplot_gene-', gene, '.pdf')))
  print(p)
  dev.off()
}

#' Get gene annotation from BioMart
#'
#' Gets annotation for the set of genes desired
#' @param GeneList Gene IDs entries. ENSMBL genes only.
#' @param datasets By default: hsapiens_gene_ensembl
#' @export
get_gene_annotation <- function(Genelist, species="hsapiens_gene_ensembl") {
  ##get gene symbols from biomart - watch out for suffixes effect on annotation retrieval!!!
  
  library(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = species)
  
  ## If we filter by hgnc_symbol it gets duplicates and it is not uniq
  ## ENSEMBL id is unique
  
  # Several entries for UNIPROT ids
  GeneSymbolsTable_full <-getBM(attributes = c('ensembl_gene_id_version','ensembl_gene_id', 
                                               'hgnc_symbol', 'description', 'gene_biotype', 
                                               'uniprot_gn_id', 'uniprotswissprot'),
                                filters = 'ensembl_gene_id', 
                                values = Genelist,
                                mart = mart)
  
  ## we will be using uniprot_swissprot
  #head(GeneSymbolsTable_full)
  #length(unique(GeneSymbolsTable_full$uniprot_gn))
  #length(unique(GeneSymbolsTable_full$uniprotswissprot))
  
  ## discard missing swisprot ids
  #library(tidyr)
  GeneSymbolsTable_swissprot_filter <- GeneSymbolsTable_full[!(is.na(GeneSymbolsTable_full$uniprotswissprot) | GeneSymbolsTable_full$uniprotswissprot==""), ]
  #head(GeneSymbolsTable_swissprot_filter)
  #dim(GeneSymbolsTable_full)
  #dim(GeneSymbolsTable_swissprot_filter)
  #length(unique(GeneSymbolsTable_swissprot_filter$uniprotswissprot))
  
  #example_df <- head(GeneSymbolsTable_swissprot_filter, 1000)
  #example_df %>% group_by(uniprotswissprot) %>% mutate(uniprot_IDs = paste0(uniprot_gn, collapse = ",")) 
  
  ## merge by , uniprot_gn for each UniProt-Swissprot identified
  df_filtered <- GeneSymbolsTable_swissprot_filter %>% 
    group_by(uniprotswissprot) %>% 
    mutate(uniprot_IDs = paste0(unique(uniprot_gn_id), collapse = ","))
  
  df_tmp <- df_filtered %>% 
    group_by(ensembl_gene_id) %>% 
    mutate(hgnc_symbol_ID = paste0(unique(hgnc_symbol), collapse = ","))
  
  df_filtered2 <- df_tmp %>% 
    group_by(ensembl_gene_id) %>% 
    mutate(uniprotswissprot_IDs = paste0(unique(uniprotswissprot), collapse = ","))
  
  dim(df_filtered2)
  head(df_filtered2)
  tail(df_filtered2)
  
  ## discard column: uniprot_gn_id, uniprot_swissprot
  df_filtered2$uniprot_gn_id <- NULL
  df_filtered2$uniprotswissprot <- NULL
  df_filtered2$hgnc_symbol <- NULL
  
  dim(df_filtered)
  head(df_filtered)
  
  ## remove duplicates
  df_filtered3 <- df_filtered2[!duplicated(df_filtered2),]
  dim(df_filtered3)
  head(df_filtered3)
  
  return(df_filtered3)
  
}


#' Get gene coordinates from BioMart
#'
#' Gets coordinates for the set of genes desired
#' @param GeneList Gene IDs entries. ENSMBL genes only.
#' @param datasets By default: hsapiens_gene_ensembl
#' @export
get_gene_coordinates <- function(Genelist, species="hsapiens_gene_ensembl") {
  ##get gene symbols from biomart - watch out for suffixes effect on annotation retrieval!!!
  
  library(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = species)
  
  ## If we filter by hgnc_symbol it gets duplicates and it is not uniq
  ## ENSEMBL id is unique
  
  # Several entries for UNIPROT ids
  Gene_coordinatesTable <-getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand'),
                                filters = 'ensembl_gene_id', 
                                values = Genelist,
                                mart = mart)
  
  return(Gene_coordinatesTable)
}

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

#' Create topGO data and enrichment analysis
#' 
#' Creates Gene Ontology data using topGO package and enrichment analysis
#' @param list_genes List of DE gene (EntrezID).
#' @param list_all_genes List of all available genes (EntrezID).
#' @param gene2go ViSEAGO annotation object for the package annotation of interest.
#' @param ont_given Ontology to analyse: Biological process (B); Mollecular function (M) or Celullar component (C) 
#' @param nodeSize_given Size of the node used to prune the GO hierarchy from the terms which have less than X annotated genes
#' @param cutoff_given Pvalue cutoff given
#' @export
get_topGO_data <- function(list_genes, list_all_genes, gene2GO, ont_given="B", nodeSize_given=5, cutoff_given=0.05, statistic_given="fisher") {
  
  library(ViSEAGO)
  library(topGO)
  
  ## create topGO data 
  topGO_data <- ViSEAGO::create_topGOdata(geneSel = list_genes, allGenes = list_all_genes,  
                                          gene2GO = gene2GO,  ont=ont_given, nodeSize = nodeSize_given)
  
  ## topGO_data: contains all gene identifiers and their scores, GO annotations, GO hierarchical structure
  ##             and additional information require to perform enrichment analysis
  ##
  ## **NOTE: argument nodeSize: is used to prune the GO hierarchy from the terms which have less than X annotated genes
  ## 
  ## **NOTE: you can create topGO data directly from topGO but here I use ViSEAGO to create it and also to retrieve gene2GO object from Bioconductor
  ##         annotation package using 
  ##    
  ##         Bioconductor <- ViSEAGO::Bioconductor2GO()
  ##         myGene2GO <- ViSEAGO::annotate("org.Mm.eg.db", Bioconductor)
  ##    
  
  ## Enrichment analysis
  ##
  ## Two types of statistics available: 
  ##  fisher: Fisher exact test: based on gene counts
  ##  ks: Kolmogorov-Smirnov like test whihc computes enrichment based on gene scores
  
  ## algorithm classic:
  ##    Tests over-representation of GO terms within the group of differentially expressed genes. Each GO category is tested independently
  classic_data <- topGO::runTest( topGO_data, algorithm = "classic", statistic = statistic_given, cutOff=cutoff_given)
  
  ## algorithm elim:
  ##    It was design to be more conservative then the classic method and therefore one expects the p-values returned to be lower.
  elim_data <- topGO::runTest( topGO_data, algorithm = "elim", statistic = statistic_given, cutOff=cutoff_given)
  
  ## return data
  list2return <- list(
    "topGO" = topGO_data,
    "elim" = elim_data,
    "classic" = classic_data
  )
  
  return(list2return)
}

#' Get GSEA Datasets 
#' 
#' Gets GSEA datasets and prepares data for FGSEA analyasis
#' @param species_given Species to retrieve data. Homo sapiens by default
#' @export
get_GSEA_datasets <- function(species_given="Homo sapiens"){

  library(msigdbr)
  
  ## ATTENTION: no specific gene set Hallmark for Mus musculus although stated in the example from CRAN
  ## Please be aware that the homologs were computationally predicted for distinct genes. The full pathways may not be well conserved across species.
  
  ## download all datasets
  
  #all_gene_sets <- msigdbr(species = "Mus musculus")
  #unique(all_gene_sets$gs_cat) ## Categories
  
  #  H: Hallmark gene sets: well defined biological states or processes
  # C1: Positional gene sets: for each human chromosome
  # C2: Curated gene sets: from online pathway databases, publications, etc
  # C3: Regulatory target gene sets: gene target predictions for miRNA seed sequences and predicted TF binding sites
  # C4: Computational gene sets: defined by mining large collections of cancer-oriented microarray data
  # C5: Ontology gene sets: genes annotated by the same ontology term
  # C6: Oncogenic sigature gene sets: defined from microarray gene expression from cancer patients
  # C7: Immunological signature gene sets: represent cell states and perturbation in immune system
  # C8: Cell type signature gene sets: curated from cluster markers identified in single-cell studies
  
  ## hallmark
  print("+ Download Hallmark gene sets...")
  hallmark_gene_sets <- msigdbr(species =species_given, category = "H") 
  
  ## immune
  print("+ Download Immune gene sets...")
  immune_gene_sets <- msigdbr(species = species_given, category = "C7") ## immune
  immunesigdb_gene_sets <- subset(immune_gene_sets, gs_subcat=="IMMUNESIGDB")
  
  ## ontology
  print("+ Download Ontology gene sets...")
  ontology_gene_sets <- msigdbr(species = species_given, category = "C5") ## ontology
  GO_BP_ontology_gene_sets <- subset(ontology_gene_sets,gs_subcat=="GO:BP")
  GO_CC_ontology_gene_sets <- subset(ontology_gene_sets,gs_subcat=="GO:CC")
  GO_MF_ontology_gene_sets <- subset(ontology_gene_sets,gs_subcat=="GO:MF")
  
  ## curated pathway databases
  print("+ Download Pathway gene sets...")
  pathway_gene_sets <- msigdbr(species = species_given, category = "C2") ## curated gene sets
  
  ## regulatory databases
  print("+ Download Regulatory gene sets...")
  regulatory_gene_sets <- msigdbr(species = species_given, category = "C3") ## regulatory
  
  ## ---------------------
  ## create list for FGSEA
  ## ---------------------
  h_msigdbr_list = split(x=hallmark_gene_sets$gene_symbol, f=hallmark_gene_sets$gs_name)
  
  i_msigdbr_list = split(x=immunesigdb_gene_sets$gene_symbol, f=immunesigdb_gene_sets$gs_name)
  
  GO_BP_msigdbr_list = split(x=GO_BP_ontology_gene_sets$gene_symbol, f=GO_BP_ontology_gene_sets$gs_name)
  GO_CC_msigdbr_list = split(x=GO_CC_ontology_gene_sets$gene_symbol, f=GO_CC_ontology_gene_sets$gs_name)
  GO_MF_msigdbr_list = split(x=GO_MF_ontology_gene_sets$gene_symbol, f=GO_MF_ontology_gene_sets$gs_name)
  
  pathway_msigdbr_list = split(x=pathway_gene_sets$gene_symbol, f=pathway_gene_sets$gs_name)
  
  regulatory_msigdbr_list = split(x=regulatory_gene_sets$gene_symbol, f=regulatory_gene_sets$gs_name)
  
  ## save as excel information from msigdb
  #hallmark_gene_sets_info <- as.data.frame(unique(hallmark_gene_sets[,c(3,15)]))
  #immunesigdb_gene_sets_info <- as.data.frame(unique(immunesigdb_gene_sets[,c(1:3,10:15)]))
  #ontology_gene_sets_info <- as.data.frame(unique(ontology_gene_sets[,c(2,3,13:15)]))
  
  list2return = list(
    "hallmark" = h_msigdbr_list,
    "immune" = i_msigdbr_list,
    "GO_BP" = GO_BP_msigdbr_list,
    "GO_CC" = GO_CC_msigdbr_list,
    "GO_MF" = GO_MF_msigdbr_list,
    "pathway" = pathway_msigdbr_list,
    "regulatory" = regulatory_msigdbr_list
  )
  
  return(list2return)
}

#' Rank table by given variable
#' 
#' Ranks a given table by the given column name. It requires a GENE_SYMBOL column to include gene IDs.
#' @param table_data Dataframe containing information
#' @param option_given Column name to retrieve data.
#' @export
rank_list_by <- function(table_data, option_given="logFC") {
  
  table_data <- table_data[!table_data$GENE_SYMBOL=="",]
  gene_list <- table_data[[option_given]]
  names(gene_list) = table_data$GENE_SYMBOL
  
  gene_list <- gene_list[order(gene_list)]
  
  return(gene_list)
}

#' FGSEA enrichment analysis
#' 
#' Creates a gene set enrichment analysis using fgsea a given table by the given column name. It requires a GENE_SYMBOL column to include gene IDs.
#' @param gene_list_provided Named list ranked by either pvalue, logFC, padj
#' @param myGeneSet Gene set data as a list of lists
#' @param title_given Title to include in the ggplot generated
#' @param nproc_given Number of threads to use
#' @export
FGSEA_GSEA <- function(gene_list_provided, myGeneSet, title_given="example", nproc_given=2, minSize=10, maxSize=600, nPermSimple = 10000) {
  
  library(fgsea)
  fgRes <- fgsea::fgseaMultilevel(pathways = myGeneSet, minSize=minSize, maxSize=maxSize, nPermSimple = nPermSimple,
                                  stats = gene_list_provided, nproc = nproc_given) %>% as.data.frame()
  #minSize=1,
  #minSize=15, maxSize=600, nperm=10000) %>% as.data.frame()
  
  fgRes1 <- fgRes[fgRes$padj<0.25,]
  data.table::setorder(fgRes1, padj)
  p <- plot_GSEA(fgRes1, title_given)
  
  returnList = list(
    fgRes = fgRes,
    p = p
  )
  
  return(returnList)
}

#' Plog GSEA
#' 
#' Plots GSEA results using ggplot and enrichment score provided
#' @param fgRes Dataframe containing information
#' @param title_given Title to include in the ggplot generated
#' @export
plot_GSEA <- function(fgRes, title_given) {
  
  fgRes$Enrichment = ifelse(fgRes$ES > 0, "Up-regulated", "Down-regulated")
  
  filtRes = rbind(head(fgRes, n = 10),tail(fgRes, n = 10 ))
  #filtRes = fgRes
  
  library(ggplot2)
  g = ggplot(filtRes, aes(reorder(pathway, ES), ES)) +
    geom_segment( aes(reorder(pathway, ES), xend=pathway, y=0, yend=ES)) +
    geom_point( size=5, aes( fill = Enrichment),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                                 "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=title_given) + 
    theme_minimal()
  
  return(g)
}

#' Save plot in pdf
#' 
#' Save a plot in a pdf in the folder provided.
#' @param folder_path Absolute path to save pdf
#' @param name_file PDF file name, do not include extension
#' @param plot_given Plot object
#' @export
save_pdf <- function(folder_path, name_file, plot_given) {
  pdf(file.path(folder_path, paste0(name_file, ".pdf")), paper = "A4r", width = 35, height = 12)
  print(plot_given)
  dev.off()
}
