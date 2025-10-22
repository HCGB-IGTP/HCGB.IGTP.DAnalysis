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
#' @param min_count integer for minimum count cutoff (Default=10). 
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


#' Discard low counts in dataframe
#'
#' This functions discard low counts in a dataframe with any non-numeric column
#' @param df_given Dataframe to pase
#' @param min_count Minimum count to use as cutoff (Default: 10).
#' @export
discard_lowCounts_df = function(df_given, min_count=10) {
  df_given <- df_given %>% mutate(total=rowSums(select_if(., is.numeric)))
  keep <- df_given['total'] >= min_count
  df_given['total'] <- NULL
  return(df_given[keep,])
}

#' Create DE gene plots for several metadata variables
#'
#' @param outdir_path Output dirname
#' @param metadata_given Samplesheet dataframe
#' @param list_of_cols Columns to use for plotting genes
#' @param sign.data Signficant dataframe to use
#' @param comp_name Comparison of interest
#' @param numerator Group 1, used as numerator
#' @param denominator Group 2, used as denominator or reference in the comparison
#' @param gene.annot.df Dataframe with annotation, either for miRNA or mRNA containing parent and variant or hgnc_symbol, respectively. Default: NULL
#' @param max_cutoff_to_plot Maximun number of genes to plot. Default: 50
#' @param data_type Either mRNA, miRNA or rRNA_16S. Default: mRNA
#' @param outdir_name Name of the folder to create. Default: boxplot_DE
#' @param input_type Choose either DESeq2 or limma to identify the format of the results.
#' @export
create_DE_plots <- function(outdir_path, metadata_given, list_of_cols, sign.data, 
                            comp_name, numerator, denominator,
                            gene.annot.df=NULL, max_cutoff_to_plot=50, data_type="mRNA", n_cores=2, input_type="DESeq2") {
  
  print("+ Create dir")
  boxplot_DE <- outdir_path
  dir.create(boxplot_DE)
  print(boxplot_DE)
  
  ## create dataframe to plot
  print("+ Create dataframe to plot")
  
  samples_here <- rownames(metadata_given)[rownames(metadata_given) %in% colnames(sign.data)]
  
  DE_plots.df <- data.frame(row.names = samples_here, 
                            metadata_given[samples_here, list_of_cols], t(sign.data[, samples_here]))
  print(head(DE_plots.df))
  
  ##
  print("+ Iterate for each gene and create plots")
  DE_plots = list()
  
  library(doParallel)
  # Register cluster
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  
  if (max_cutoff_to_plot > length(rownames(sign.data))) {
    max_cutoff_to_plot = length(rownames(sign.data))
  }
  
  results <- foreach(gene_n = 1:max_cutoff_to_plot) %dopar% {
    gene_given = rownames(sign.data)[gene_n]
    g <- gsub("-", "\\.", gene_given)
    g <- gsub("&", "\\.", g)
    g <- gsub(":", "\\.", g)
    g <- gsub("\\+", "\\.", g)
    print(g)
    
    ## If annotation provided, add description and subtitle
    if (!is.null(gene.annot.df)) {
      gene_annot.df <- gene.annot.df[gene_given, ]
      if (data_type == "mRNA") {
        gene_name_file = paste0(gene_given, "_", gene_annot.df$hgnc_symbol)
        gene_name = paste0(gene_given, "_", gene_annot.df$hgnc_symbol)
      }
      else if (data_type == "miRNA") {
        gene_name = paste0(gene_annot.df$parent, "_", 
                           gene_annot.df$variant, "_", gene_given)
        gene_name_file = paste0(gene_annot.df$parent, 
                                "_", gene_given)
      } else if (data_type == "rRNA_16S"){
        gene_name = paste0(gene_annot.df$Phylum, "_", 
                           gene_annot.df$Genus, "_", gene_given)
        gene_name_file = paste0(gene_annot.df$Genus, 
                                "_", gene_given)
      }
      ##
      print(paste0("Annotation added: ", gene_name_file))
    }
    else {
      gene_name = g
      gene_name_file = g
    }
    
    DE_plots[[gene_name_file]] = list()
    for (i in colnames(metadata_given[, list_of_cols])) {
      g <- gsub("-", "\\.", g)
      print(paste0("Variable: ", i))
      if (is.numeric(metadata_given[, i])) {
        print(paste0("Numeric: ", i))
        p2.tmp <- HCGB.IGTP.DAnalysis::ggscatter_plotRegression(data_all_given = DE_plots.df, 
                                                                x.given = g, y.given = i, title_string = i)
        p2 <- p2.tmp$plot + ggtitle(paste0("Variable: ", i), subtitle = gene_name)
      }
      else {
        p2 <- HCGB.IGTP.DAnalysis::ggboxplot_scatter(data_all_given = DE_plots.df, 
                                                     colName = i, y.coord = g)
        
        
        if (i == comp_name) {
          ## If the plot is for the comparison of interest, 
          
          if (input_type=="DESeq2") {
            
            ## add in blue the pvalue adjusted obtained by DESeq2
            stat.test <- tibble::tribble(~group1, ~group2, 
                                         ~p.adj, numerator, denominator,
                                         sign.data[gene_given, "padj"])
            p2 <- p2 + stat_pvalue_manual(stat.test, 
                                          y.position = max(p2$data[,g]) * 0.85, 
                                          step.increase = 1, color = "blue", 
                                          label = "p.adj")
            
            p2 <- p2 + labs(title = paste0("Variable: ", i),
                            subtitle = gene_name,
                            caption = "[**In blue, pvalue adjusted by DESeq2]")  
            
          } else if (input_type=="limma") {
            ## add in blue the pvalue adjusted obtained by limma
            stat.test <- tibble::tribble(~group1, ~group2, 
                                         ~p.adj, numerator, denominator,
                                         sign.data[gene_given, "adj.P.Val"])
            p2 <- p2 + stat_pvalue_manual(stat.test, 
                                          y.position = max(p2$data[,g]) * 0.85, 
                                          step.increase = 1, color = "blue", 
                                          label = "p.adj")
            
            p2 <- p2 + labs(title = paste0("Variable: ", i),
                            subtitle = gene_name,
                            caption = "[**In blue, pvalue adjusted by limma]")  
            
          }
          
          
        } else {
          p2 <- p2 + ggtitle(paste0("Variable: ", i), subtitle = gene_name)  
        }
        
      }
      print(p2)
      DE_plots[[gene_name_file]][[i]] = p2
    }
    
    HCGB.IGTP.DAnalysis::save_multi_pdf(folder_path = boxplot_DE, name_file = gene_name_file, list_plots = DE_plots[[gene_name_file]])
  }
  
  # Don't fotget to stop the cluster
  stopCluster(cl = cluster)
  
  return(DE_plots.df)
  
}


#' Filter significant hits from differential expression analysis
#' 
#' When running DESeq2 you usually required significant hits.
#' @param dataF Dataframe with either normalized values and DESEQ2 values or only DESEQ2.
#' @param sign_value Pvalue adjusted cutoff: Default=0.05
#' @param LFC Log Fold Change cutoff: Default: 0.26
#' @param input_type Choose either DESeq2 or limma to identify the format of the results.
#' @export
filter_signficant_hits <- function(dataF, sign_value = 0.05, LFC=0.26, input_type="DESeq2") {
  
  if (input_type=="DESeq2") {
    #log2FoldChange
    #padj
    dataFilt <- subset(dataF, abs(log2FoldChange)>LFC & padj<sign_value)
    dataFilt <- dataFilt[order(dataFilt$padj),]  
  
  } else if (input_type=="limma") {
    #log2FoldChange
    #padj
    dataFilt <- subset(dataF, abs(log2FoldChange)>LFC & padj<sign_value)
    dataFilt <- dataFilt[order(dataFilt$padj),]
  }

  
  return(dataFilt)
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


