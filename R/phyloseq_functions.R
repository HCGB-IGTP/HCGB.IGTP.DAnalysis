#' Ordination plots
#'
#' @param dist_method phyloseq::distance method to use (Available in physeq::distanceMethodList())
#' @param phyloseq_obj Phyloseq object to plot
#' @param color_given Metadata column to use for color
#' @param shape_given Metadata column to use for shape
#' @param ordination_method Ordination method. Options: c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"). Default: MDS
#'
#' @export
distance_ord_plot <- function(dist_method, phyloseq_obj, color_given, shape_given, ordination_method="MDS") {
  
  library(phyloseq)
  library(ggplot2)
  
  print("Method:")
  print(dist_method)
  
  if (dist_method %in% unlist(phyloseq::distanceMethodList) ) {
    
  } else {
    print("ERROR: Not avaialble")
    print("Find options in: phyloseq::distanceMethodList")
    return()
  }
  
  if (dist_method == "bray") {
    # Transform data to proportions as appropriate for Bray-Curtis distances
    phyloseq_obj <- transform_sample_counts(phyloseq_obj, function(otu) otu/sum(otu))
  }
  
  
  # Calculate distance matrix
  iDist <- phyloseq::distance(physeq = phyloseq_obj, method=dist_method)
  # Calculate ordination
  iORD  <- phyloseq::ordinate(phyloseq_obj, method = ordination_method, distance=iDist)
  
  ## Make plot
  # Create plot, store as temp variable, p
  p <- phyloseq::plot_ordination(phyloseq_obj, ordination = iORD, color=color_given, shape=shape_given)
  
  # Add title to each plot
  p <- p + ggtitle(paste(ordination_method, " using distance method ", dist_method, sep=""), 
                   subtitle = paste0("Color: ", color_given, "; Shape: ", shape_given)) + 
    geom_point(size = 3) + scale_shape_manual(values = c(3,17)) + theme_classic()
  
  # return
  p
}

#' Get most abundant taxa
#'
#' @param physeq_obj Phyloseq object to subset
#' @param n Number of top taxa to retrieve
#'
#' @returns
#' @export
#'
#' @examples
get_abundant_taxa <- function(physeq_obj, n=10) {

  library(phyloseq)
  
  # Get the most abundant N taxa
  TopNOTUs <- names(sort(taxa_sums(physeq_obj), TRUE)[1:n]) 
  print(paste0("+ Prune taxa. n=", n ))
  
  physeq.top   <- prune_taxa(TopNOTUs, physeq_obj)
  
  ## transform sample counts to get relative abundance
  physeq.top.prop <- transform_sample_counts(physeq.top, function(OTU) OTU/sum(OTU))
  
  return(list("top"=physeq.top,
              "top.prop"=physeq.top.prop))
  
}

#' Get richness & diversity statistics
#'
#' @param physeq.obj Phyloseq object to use
#' @param results.fold Results folder to store data in csv and xlsx format
#'
#' @returns dataframe
#' @export
get_diversity_stats <- function(physeq.obj, results.fold) {
  
  alpha_measures_diversity = c("Shannon", "Simpson", "InvSimpson", "Fisher")
  alpha_measures_richness = c("Observed", "ACE", "Chao1")
  
  both_measures <- c(alpha_measures_diversity, alpha_measures_richness)
  
  alphaDiversity_diversity <-  estimate_richness(physeq.obj, measures = alpha_measures_diversity)
  alphaDiversity_richness <-  estimate_richness(physeq.obj, measures = alpha_measures_richness)
  
  ## save values to csv
  write.csv(alphaDiversity_diversity, file.path(results.fold, "alphaDiversity_diversity.csv"))
  write.csv(alphaDiversity_richness, file.path(results.fold, "alphaDiversity_richness.csv"))
  head(alphaDiversity_diversity) 
  head(alphaDiversity_richness) 
  
  ## merged all data
  merged_data <- merge(physeq.obj@sam_data, 
                       cbind(alphaDiversity_diversity, alphaDiversity_richness), 
                       by="row.names")
  head(merged_data)
  rownames(merged_data) <- merged_data$Row.names
  merged_data$Row.names <- NULL
  
  HCGB.IGTP.DAnalysis::save_woorkbook_openxlsx(file_path = file.path(results.fold, "alphaDiversity_richness.xlsx"), 
                                               list_df = list("alphaDiversity"=merged_data))
  
  return(merged_data)
  
}


#' Get abundances in a long format and save file in folder
#' 
#' @param physeq.obj Phyloseq object to use
#' @param results.fold Results folder to store data in csv and xlsx format
#' @param tag_name Tag to add to the name of the files generated
#' @export
get_abundance.long <- function(physeq.obj, results.fold, tag_name="data") {
  # ## get abundances melted
  long_abundances <- reshape2::melt(otu_table(physeq.obj), measure.vars=rownames(otu_table(physeq.obj)), 
                           variable.name="sample", 
                           value.name="proportion")
  head(long_abundances)
  write.csv(long_abundances, file = file.path(results.fold, paste0(tag_name, "_abundance.csv")))
  
  return(long_abundances)
}


#' Get prevalence 
#'
#' @param physeq.obj Phyloseq object to use
#' @param results.fold Results folder to store data in csv and xlsx format
#' @param tag_name Tag to add to the name of the files generated
#' @param prevThres Proportion of samples as threshold
#'
#' @export
get_prevalence <- function(physeq.obj, results.fold, tag_name="data", prevThres=0.05) {
  
  library(phyloseq)
  prevdf = apply(X = otu_table(physeq.obj),
                 MARGIN = ifelse(taxa_are_rows(physeq.obj), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  
  # Le agregamos la taxonomía
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq.obj),
                      tax_table(physeq.obj))
  
  dfprev <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq.obj, "Phylum"))
  prevalence_plot <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq.obj),color=Phylum)) +
    # Agregamos una línea para nuestro umbral
    geom_hline(yintercept = prevThres, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  
  print(prevalence_plot)
  
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = results.fold, 
                                name_file = paste0(tag_name, "prevalence_plot"), 
                                plot_given = prevalence_plot)
  
  prevalenceThreshold = prevThres * nsamples(physeq.obj)
  
  keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
  physeq.filt = prune_taxa(keepTaxa, physeq.obj)
  
  return(list("plot"=prevalence_plot, 
              "physeq.obj"=physeq.filt))
  
  
}



#' Geometric mean
#'
#' Calculate geometric means prior to estimate size factors
#' @param x Matrix
#' @param na.rm Remove NAs
#'
#' @export
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#' Call DESeq2 from phyloseq object
#'
#' @param physeq.obj Phyloseq object to use
#' @param list_of_cols List of columns of interest to subset from metadata
#' @param data_type Either mRNA, miRNA or rRNA_16S. Default: rRNA_16S
#' @param tax_levels Taxonomic levels to use. Default: c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' @param sample_sheet_given Samplesheet containing metadata information
#' @param OUTPUT_Data_dir_given Absolute path to store results
#' @param int_threads Number of threads to use
#' @param formula_given Design formula to use
#' @param coef_n Number of the coefficient of results to test (if desired)
#' @param early_return Whether to return exploratory results early or not
#' @param comp_ID Tag name to include for each comparison
#' @param cutoff.given add an option to include cutoff when removing Zeros
#' @param sign_value.given Adjusted pvalue cutoff. Default=0.05, 
#' @param LFC.given Log Fold change cutoff. Default=log2(1.2), 
#' @param localFit Use a fitType=local for mean dispersion fit in DESeq2
#' @param forceResults Boolean to force re-run analysis if already generated in the folder provided
#' @param gene.annot Dataframe containing gene annotation (Default: NULL)
#' @param shrinkage.given LFC shrinkage estimator provided. Available: apeglm, ashr or normal
#' @param min_cutoff_to_plot Minimun number of genes significant to continue analysis. Default=3
#' @param max_cutoff_to_plot Number of genes significant to plot as candidates analysis. Default=50
#' @export
analysis_DESeq_phyloseq <- function(OUTPUT_Data_dir_given, physeq.obj, sample_sheet_given, 
                                    list_of_cols, formula_given, int_threads=2,
                                    sign_value.given = 0.05, LFC.given = log2(1.2),
                                    coef_n=NA, early_return=FALSE, comp_ID=NULL, cutoff.given=0.9, 
                                    localFit=FALSE, forceResults=FALSE, min_cutoff_to_plot=3, max_cutoff_to_plot=50,
                                    gene.annot=NULL, data_type = "rRNA_16S",
                                    shrinkage.given="apeglm", tax_levels = c("Domain", "Phylum", "Class", 
                                                                             "Order", "Family", "Genus", 
                                                                             "Species")) {
  
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
  library(phyloseq)
  #--------------------------
  
  dir.create(OUTPUT_Data_dir_given, showWarnings = FALSE)
  
  ## Set parallel threads
  print (paste0("Set Multicore: ", int_threads))
  register(MulticoreParam(int_threads))
  
  #############
  ## Design ###
  #############
  ddsFullCountTable = phyloseq_to_deseq2(physeq.obj, design = as.formula(formula_given))
  geoMeans = apply(counts(ddsFullCountTable), 1, gm_mean)
  ddsFullCountTable = estimateSizeFactors(ddsFullCountTable, geoMeans = geoMeans)
  
  ############ 
  # Norm.counts
  ############
  # Save normalized values
  # Normalized values
  
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
  
  dds_object1 = HCGB.IGTP.DAnalysis::discard_lowCounts(dds_object = dds_object)
  
  normValues <- counts(dds_object1, normalized=T)
  print(head(normValues))
  
  
  
  if (!is.null(gene.annot)) {
    
    print("Remove some characters that might create errors in the taxonomy table...")
    
    ## remove some characters that might create errors
    ## Not only in genes, in all other columns
    gene.annot <- gene.annot %>% dplyr::mutate(across(tax_levels, stringr::str_remove_all, "/"))
    print(head(gene.annot))
    print("...")
    
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
                                 sign_value.given = sign_value.given, LFC.given = LFC.given, 
                                 min_cutoff_to_plot=min_cutoff_to_plot, max_cutoff_to_plot=max_cutoff_to_plot,
                                 gene.annot=gene.annot, data_type = data_type,
                                 int_threads = int_threads, coef_n = coef_n, forceResults=forceResults, shrinkage.given=shrinkage.given)
  #############
  
  #############
  ## Init data to return
  #############
  data2return <- list(
    "dds_obj" = dds_object,
    "exploratory_plots" = exploratory_plots_returned,
    "resultsNames" = resultsNames(dds_object),
    "formula" = formula_given,
    "results" = results_list
  )
  
  return(data2return)
}



