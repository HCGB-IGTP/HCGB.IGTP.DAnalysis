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
dreduction_function <- function(rnb.set_given, vars_dreduction, methods_dreduction, pdf_file) {
  pdf(pdf_file)
  for (i in vars_dreduction) {
    print (paste0("Variable: ", i))
    
    for (f in methods_dreduction) {
      print (paste0("Reduction method: ", f))
      ## generate plot Dreduction
      p<- rnb.plot.dreduction(rnb.set_given, dimensions = c(1, 2), 
                              plot.type = f,
                              point.types = "Sample_Name", point.colors=i)  
      
      p2 <- p + ggtitle(paste0("Variable: ", i,". Method: ", f))
      print(p2)
    }
  }
  dev.off()
}
###
beta_values_plots <- function(rnb.set_given, variablesOfInterest, pdf_name) {
  ## calculate methylation values
  print("Calculate methylation values")
  mm_table <- meth(rnb.set_given)  
  
  ## calculate distribution and generate plots
  print ("Generate plots")
  pdf(pdf_name)
  for (i in variablesOfInterest) {
    
    print (paste0("Variable: ",i))
    p <- rnb.plot.betadistribution.sampleGroups(beta.matrix = mm_table, 
                                                sample.group.inds = rnb.sample.groups(rnb.set_given)[[i]], 
                                                annotation = i)
    p2 <- p + ggtitle(paste0("Beta Distribution group by :", i))
    print(p2)
  }
  dev.off()
}
### 
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
###
my_volcanoplot_region <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topleft", textcx=1, point_cex = 0.3, ...) {
  
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
###
my_volcanoplot_sites <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topleft", textcx=1, point_cex = 0.3, ...) {
  
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
