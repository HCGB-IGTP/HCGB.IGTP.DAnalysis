library(RnBeads)

#' Saves RnBead object
#'
#' This functions saves a given RnBeads object in a file and dir provided
#' @param RnBeads_object RnBeads object to save
#' @param file2write File name to use
#' @param dir2save Absolute path to folder to save item
#' @export
save_RnBeads_object <- function(RnBeads_object, file2write, dir2save) {
  ##
  ## Save `RnBSet` objects to disk.
  ##
  ## creates a folder
  rnbsSaveDir <- file.path(dir2save, file2write)
  save.rnb.set(RnBeads_object, rnbsSaveDir, archive=FALSE)
  
  ## returns path to folder
  return(rnbsSaveDir)
}
##
#' Loads RnBead object
#'
#' This functions loads a previously saved RnBeads object in a file and dir.
#' 
#' @param dir4loading Absolute path to folder to load RnBeads item
#' @export
load_RnBeads_object <- function(dir4loading) {
  ##
  ## Load saved `RnBSet` objects from disk.
  ##
  rnbs.reloaded <- load.rnb.set(dir4loading)
  ## returns object loaded
  return(rnbs.reloaded)
}

##
#' Check for covariates using RnBead object
#'
#' This functions checks for covariates at cpg and promoters
#' 
#' @param var2check Suspected variable to test for Covariate
#' @param rnb.set_given RnBeads rnb.set object
#' @param plot_dir Folder to save plot results
#' @export
check_covariates <- function(var2check, rnb.set_given, plot_dir) {
  reg.types <- c("cpgislands", "promoters")
  
  # compute differential methylation tables: unadjusted and SVA adjusted
  print ("Computing differential methylation tables: unadjusted ")
  diffmeth.base <- rnb.execute.computeDiffMeth(
    x=rnb.set_given, pheno.cols=var2check, region.types=reg.types,
    adjust.sva=FALSE
  )
  
  print ("...adjusting for SVA variables")
  diffmeth.sva <- rnb.execute.computeDiffMeth(
    x=rnb.set_given, pheno.cols=var2check, region.types=reg.types,
    adjust.sva=TRUE, pheno.cols.adjust.sva=var2check
  )
  
  ## print components
  list_cmp <- as.character(get.comparisons(diffmeth.base))
  print (list_cmp)
  #get.comparisons(diffmeth.sva)
  
  ##
  pdf(paste0(plot_dir, "/Covariate_plot_", var2check, ".pdf"))
  for(i in list_cmp) {
    my_comp <- i
    print (paste0("Analyzing: ", my_comp))
    print ("Get information")
    dm.tab.base <- get.table(diffmeth.base, comparison=my_comp, "sites", return.data.frame=TRUE)
    dm.tab.sva <- get.table(diffmeth.sva,  comparison=my_comp, "sites", return.data.frame=TRUE)
    
    ## remove NAs
    print ("Remove NAs")
    dm.tab.base <- na.omit(dm.tab.base)
    dm.tab.sva <- na.omit(dm.tab.sva)
    
    # compute quantiles of -log10 p-values and prepare a data.frame for plotting
    print ("Compute quantiles of -log10 p-values")
    p.val.perc.base <- quantile(-log10(dm.tab.base$diffmeth.p.val),probs = seq(0, 1, 0.01))
    p.val.perc.sva <- quantile(-log10(dm.tab.sva$diffmeth.p.val),probs = seq(0, 1, 0.01))
    
    df <- data.frame(
      neg.log10.p.val.base=p.val.perc.base,
      neg.log10.p.val.sva=p.val.perc.sva,
      quantile=seq(0, 1, 0.01)
    )
    # plot using ggplot2
    print("Plotting")
    library(ggplot2)
    p<- ggplot(df,
               aes(x=neg.log10.p.val.base,
                   y=neg.log10.p.val.sva,
                   color=quantile)) +
      geom_point() + 
      geom_abline(intercept=0, slope=1) + 
      xlim(0, 5) + 
      ylim(0, 5) +
      ggtitle(my_comp)
    
    print(p)
    
  }
  dev.off()
}

##
#' Creates dimension reduction plots using RnBead object
#'
#' This functions creates dimension reduction (PCA & MDS) for RnBead co
#' 
#' @param rnb.set_given RnBeads rnb.set object
#' @param vars_dreduction Suspected variables to include in plot
#' @param methods_dreduction Methods for dimension reduction
#' @param plot_dir Folder to save plot results
#' @param point_types_column Column in annotation dataframe to use for naming samples.
#' @export
dreduction_function <- function(rnb.set_given, vars_dreduction, methods_dreduction, 
                                point_types_column, plot_dir) {
  
  dir.create(plot_dir)
  
  for (i in vars_dreduction) {
    print (paste0("Variable: ", i))
    
    for (f in methods_dreduction) {
      pdf_file <- paste0(i, "-", f, ".pdf")
      print (paste0("Reduction method: ", f))
      ## generate plot Dreduction
      p<- rnb.plot.dreduction(rnb.set_given, dimensions = c(1, 2), 
                              plot.type = f, 
                              point.types = point_types_column, 
                              point.colors=i)  
      
      p2 <- p + ggtitle(paste0("Variable: ", i,". Method: ", f))
      save_pdf(folder_path = plot_dir, name_file = pdf_file, plot_given = p2)
    }
  }
}

##
#' BetaDistribution plots
#'
#' This functions creates beta distribution plots for a set of variables of interest
#' 
#' @param rnb.set_given RnBeads rnb.set object
#' @param variablesOfInterest Suspected variables to include in plot
#' @param pdf_file_folder Folder to save plot results
#' @export
beta_values_plots <- function(rnb.set_given, variablesOfInterest, pdf_file_folder) {
  ## calculate methylation values
  
  dir.create(pdf_file_folder)
  
  print("Calculate methylation values")
  mm_table <- meth(rnb.set_given)  
  
  ## calculate distribution and generate plots
  print ("Generate plots")
  for (i in variablesOfInterest) {
    
    print (paste0("Variable: ",i))
    p <- rnb.plot.betadistribution.sampleGroups(beta.matrix = mm_table, 
                                                sample.group.inds = rnb.sample.groups(rnb.set_given)[[i]], 
                                                annotation = i)
    p2 <- p + ggtitle(paste0("Beta Distribution group by :", i))
    
    save_pdf(folder_path = pdf_file_folder, 
                                  name_file = i, plot_given = p2)
  }
}

##
#' Loads custom annotation 
#'
#' This functions loads and creates custom annotation from IGTP server
#' 
#' @export
custom_annotation <- function() {
  # Loading our custom enhancer annotations
  
  # The *.bed file with the annotations must be placed at the current working directory
  # Enhancers from encode
  EncodeEnhancers <- read.table('/imppc/labs/lslab/share/restored_data/imppc/labs/lslab/share/data/proc_data/IDIBELL_Met450_14_reanalysis_180417/hg19_Encode_Valid_only.bed', sep = '\t')
  colnames(EncodeEnhancers) <- c("Chromosome", "Start", "End", "name")
  rnb.set.annotation("EncodeEnhancers", EncodeEnhancers, assembly="hg19")
  
  # Enhancers from encode + FANTOM5 project
  AllEnhancers <- read.table('/imppc/labs/lslab/share/restored_data/imppc/labs/lslab/share/data/proc_data/IDIBELL_Met450_14_reanalysis_180417/hg19_FINAL_Encode_FANTOM5_enhancers.bed', sep = '\t')
  colnames(AllEnhancers) <- c("Chromosome", "Start", "End", "name")
  rnb.set.annotation("AllEnhancers", AllEnhancers, assembly="hg19")
}

##
#' VolcanoPlot region
#'
#' This functions creates volcano plot using a region of interest
#' 
#' @export
my_volcanoplot_region <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topleft", textcx=1, point_cex = 0.3, ...) {
  
  ## TODO: use Enhanced volcano
  
  with(res, plot( mean.mean.diff, -log10( comb.p.adj.fdr), pch=1, main=main, cex=point_cex, ...))
  
  with(subset(res, comb.p.adj.fdr<sigthresh ), 
       points(mean.mean.diff, -log10(comb.p.adj.fdr), pch=1, col="red", cex=point_cex, ...))
  
  with(subset(res, abs(mean.mean.diff)>lfcthresh), 
       points(mean.mean.diff, -log10(comb.p.adj.fdr), pch=1, col="orange", cex=point_cex, ...))
  
  with(subset(res, comb.p.adj.fdr<sigthresh & abs(mean.mean.diff)>lfcthresh), 
       points(mean.mean.diff, -log10(comb.p.adj.fdr), pch=1, col="green", cex=point_cex, ...))
  
  #with(subset(res, comb.p.adj.fdr<sigthresh & abs(mean.mean.diff)>lfcthresh), 
  #   textxy(mean.mean.diff, -log10(comb.p.adj.fdr), labs='', cex=0.7*textcx, ...))
  
  legend(legendpos, xjust=1, yjust=1, 
         legend=c(paste("FDR<",sigthresh,sep=""), paste("|FC|>",round(2^lfcthresh,1),sep=""), "both"), 
         pch=1, cex=.8, col=c("red","orange","green"))
}

##
#' VolcanoPlot site
#'
#' This functions creates volcano plot using a site of interest
#' 
#' @export
my_volcanoplot_sites <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topleft", textcx=1, point_cex = 0.3, ...) {
  
  ## TODO: use Enhanced volcano
  
  with(res, plot( mean.diff, -log10( diffmeth.p.adj.fdr), pch=1, main=main, cex=point_cex, ...))
  
  with(subset(res, diffmeth.p.adj.fdr<sigthresh ), 
       points(mean.diff, -log10(diffmeth.p.adj.fdr), pch=1, col="red", cex=point_cex, ...))
  
  with(subset(res, abs(mean.diff)>lfcthresh), 
       points(mean.diff, -log10(diffmeth.p.adj.fdr), pch=1, col="orange", cex=point_cex, ...))
  
  with(subset(res, diffmeth.p.adj.fdr<sigthresh & abs(mean.diff)>lfcthresh), 
       points(mean.diff, -log10(diffmeth.p.adj.fdr), pch=1, col="green", cex=point_cex, ...))
  
  #with(subset(res, diffmeth.p.adj.fdr<sigthresh & abs(mean.diff)>lfcthresh), 
  #   textxy(mean.diff, -log10(diffmeth.p.adj.fdr), labs='', cex=0.7*textcx, ...))
  
  legend(legendpos, xjust=1, yjust=1, 
         legend=c(paste("FDR<",sigthresh,sep=""), paste("|FC|>",round(2^lfcthresh,1),sep=""), "both"), 
         pch=1, cex=.8, col=c("red","orange","green"))
}

##
#' RnBeads Locus plotter caller
#'
#' This functions calls rnb.plot.locus.profile for a gene or region of interest
#' @param RnBeads_object RnBeads object
#' @param grouping_class rnb.sample.groups output grouping samples
#' @param gene_name Name of the gene/region to create pdf
#' @param comparison Name to add to create pdf
#' @param outdir Folder to store pdf
#' @param chr Chromosome ID
#' @param start Start position
#' @param end End position
#' @export
locus_plotter <- function(RnBeads_object, grouping_class, gene_name, comparison, outdir, chr, start, end) {
  
  ## plotting
  pdf_name <- file.path(outdir, paste0(gene_name, '-', comparison, '.pdf'))
  pdf(pdf_name)
  rnb.plot.locus.profile(RnBeads_object, chr, start, end, grps = grouping_class)
  dev.off()
}

##
#' RnBeads Locus plotter dataframe
#'
#' This functions calls locus_plotter (rnb.plot.locus.profile caller) for each entry in a dataframe
#' @param data_frame_given Dataframe with information from differential_methylation_data/diffMethTable_region_genes...
#' @param name_given name to use
#' @param RnBeads_object RnBeads object
#' @param sample.grouping rnb.sample.groups output grouping samples
#' @param out_folder Folder to store pdfs
#' @export
plot_genes <- function(data_frame_given, name_given, RnBead_object, sample.grouping, out_folder) {
  
  ## data_frame_given contains as columns:
  ##   "id", "Chromosome", "Start", "End", "symbol", "mean.mean.diff", "comb.p.val", "comb.p.adj.fdr"
  
  ## Iterate over a dataframe and plot methylation values  
  for (row in 1:nrow(data_frame_given)) {
    symbol_name <- data_frame_given[row, "symbol"]
    row_id <- data_frame_given[row, "id"]
    
    if (is.na(symbol_name)) {
      gene_name <- row_id
    } else {
      gene_name <- paste0(symbol_name, '-', row_id)
    }
    
    ## 
    print ("##################")
    print (gene_name)  
    print (data_frame_given[row,])
    
    locus_plotter(RnBeads_object = RnBead_object, grouping_class = sample.grouping, 
                  gene_name = gene_name, 
                  comparison = name_given, 
                  outdir = out_folder,
                  chr = data_frame_given[row, 'Chromosome'],
                  start = data_frame_given[row, 'Start'],
                  end = data_frame_given[row, 'End']
    )
    print ('')
    
  }
}

#' Load Beta Values
#'
#' This functions Reads into a dataframe the Beta values CSV file from RnBeads Track and table export module. 
#' Removes missing rows using delete_na()
#' @param betas_table_file Absolute path to CSV file
#' @param allow_missing_p Default=0.95. Allows only 5% missing data
#' @export
load_beta_values <- function(betas_table_file, allow_missing_p=0.95) {
  betas_table <- read.table(betas_table_file, sep = ',', 
                            header = TRUE, check.names = FALSE)
  
  ## remove unnecessary
  betas_table$Chromosome <- NULL
  betas_table$Strand <- NULL
  betas_table$Start <- NULL
  betas_table$End <- NULL
  
  ## rename column names
  #colnames(betas_table)[2] <- "chromosome"
  #colnames(betas_table)[3] <- "start"
  #colnames(betas_table)[4] <- "end"
  #colnames(betas_table)[5] <- "strand"
  
  ## Set rownames
  rownames(betas_table) <- betas_table$ID
  betas_table$ID <- NULL
  
  ## Allow only allow_missing_p% missing values
  betas_table <- delete_na(betas_table, round(allow_missing_p*ncol(betas_table)))
  return(betas_table)
}
