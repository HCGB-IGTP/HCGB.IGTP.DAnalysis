## original code from: ~/proc_data/20231213_EPIC_MAP/code/2.2.epic19_exploratorClinical_analysis.R

wgcna_soft_thresh <- function(input_data.mat, title2include="",
                              start_i=1, end_i=20, iter_steps=1, blockSize.given=30 ) {
  
  # The soft thresholding, is a value used to power the correlation of the genes to that threshold. 
  # The assumption on that by raising the correlation to a power will reduce the noise of the 
  # correlations in the adjacency matrix.
  # To pick up one threshold use the pickSoftThreshold function, which calculates for each power 
  # if the network resembles to a scale-free graph. The power which produce a higher similarity with a 
  # scale-free network is the one you should use.
  
  # Pick a soft threshold power near the curve of the plot and then,create the network using the blockwiseModules command. 
  
  powers = c(seq(from = start_i, to = end_i, by = iter_steps))
  
  library(WGCNA)
  
  # Call the network topology analysis function
  sft_data = pickSoftThreshold(
    data = input_data.mat,
    blockSize = blockSize.given,
    powerVector = powers,
    verbose = 5
  )
  
  library(ggplot2)
  ## scale free topology
  p1 <- ggplot(data = sft_data$fitIndices, mapping = aes(x=Power, y= -sign(slope)*SFT.R.sq)) + 
    geom_line(color="blue") + geom_point() + ggtitle(title2include) + 
    ylab("Scale Free Topology Model Fit, signed R^2") + xlab("Soft Threshold (power)") + 
    theme_classic() + geom_hline(yintercept=0.9, linetype="dashed", color = "red") + 
    ggrepel::geom_label_repel(mapping = aes(x=Power, y= -sign(slope)*SFT.R.sq, label=Power))
  
  # Mean connectivity as a function of the soft-thresholding power
  p2 <- ggplot(data = sft_data$fitIndices, mapping = aes(x=Power, y= mean.k.)) + 
    geom_line(color="blue") + geom_point() + ggtitle(title2include) + 
    ylab("Mean Connectivity") + xlab("Soft Threshold (power)") + 
    theme_classic() + geom_hline(yintercept=0.9, linetype="dashed", color = "red") 
  
  data2retuns <- list(
    "data"=sft_data,
    "plot_scale"=p1,
    "mean_connectivity"=p2
  )
}

wgcna_blockwiseModules <- function(picked_power, input_data.mat, tag_name, plot_dir) {
  # The blockwiseModule may take a while to run, 
  # since it is constructing the TOM (topological overlap matrix) and several other steps. While it runs, take a 
  # look at the blockwiseModule documentation (link to vignette) for more information on the parameters. 
  
  
  #--------
  # Turn data expression into topological overlap matrix
  #--------
  library(WGCNA)
  
  # Automatic
  temp_cor <- cor       
  cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
  netwk_res <- WGCNA::blockwiseModules(input_data.mat,             # <= Input data transposed
                                       
                                       # == Adjacency Function ==
                                       power = picked_power,                # <= power here
                                       networkType = "signed",
                                       
                                       # == Tree and Block Options ==
                                       deepSplit = 2,
                                       pamRespectsDendro = F,
                                       
                                       # detectCutHeight = 0.75,
                                       minModuleSize = 30,
                                       maxBlockSize = 4000,
                                       
                                       # == Module Adjustments ==
                                       reassignThreshold = 0,
                                       mergeCutHeight = 0.25, # 0.25
                                       
                                       # == TOM == Archive the run results in TOM file (saves time)
                                       saveTOMs = T,
                                       saveTOMFileBase = "ER",
                                       
                                       # == Output Options
                                       numericLabels = T,
                                       verbose = 3)
  
  cor <- temp_cor     # Return cor function to original namespace
  #--------
  
  # Convert labels to colors for plotting
  mergedColors = labels2colors(netwk_res$colors)
  
  #--------
  # Plot the dendrogram and the module colors underneath
  #--------
  pdf(file.path(plot_dir, paste0(tag_name, "_plotDendroAndColors.pdf")), paper = "A4r", width = 35, height = 12)
  WGCNA::plotDendroAndColors(
    netwk_res$dendrograms[[1]],
    mergedColors[netwk_res$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main= paste("Cluster dendogram - mergeCutHeight:0.25"))
  dev.off()
  
  module_df <- data.frame(
    gene_id = names(netwk_res$colors),
    colors = labels2colors(netwk_res$colors)
  )
  
  print(table(module_df$colors))
  #--------
  
  ###
  data2return <- list(
    "netwk_res" = netwk_res,
    "module_df" = module_df,
    "mergedColors" = mergedColors
  )
  
  return(data2return)
}

sample_sheet_splitter <- function(sample_sheet_given, max_cols=15) {
  max_len <- length(colnames(sample_sheet_given))
  f <- unique(cut(1:max_len,
                  breaks= round(max_len/max_cols)))
  labs <- levels(f)[f]
  lower <- round(as.numeric( sub("\\((.+),.*", "\\1", labs)))
  
  list_of_cols <- list()  
  for (i in lower) {
    letter2paste <- toupper(letters[which(lower==i)])
    list_of_cols[[letter2paste]] <- colnames(sample_sheet_given)[i:(i+max_cols)]
  }
  
  return(list_of_cols)
  
}

wgcna_eigengenes_module <- function(input_data.mat, colors2use, sample_sheet_given, plot_dir, tag_name) {
  library(WGCNA)
  
  # Get Module Eigengenes per cluster
  MEs0_res <- WGCNA::moduleEigengenes(input_data.mat, 
                                      colors2use)$eigengenes
  
  # Reorder modules so similar modules are next to each other
  MEs0_res <- orderMEs(MEs0_res)
  head(MEs0_res)
  
  ##
  print("rownames(sample_sheet_given) == rownames(MEs0_res)")
  print(rownames(sample_sheet_given) == rownames(MEs0_res))
  nGenes = ncol(input_data.mat)
  nSamples = nrow(input_data.mat)
  nums <- unlist(lapply(sample_sheet_given, is.numeric), use.names = FALSE)  
  
  moduleTraitCor_res = cor(MEs0_res, sample_sheet_given, use= "p")
  moduleTraitPvalue_res = corPvalueStudent(moduleTraitCor_res, nSamples)
  
  write.csv(moduleTraitPvalue_res, file = file.path(res_dir, paste0(tag_name, "_moduleTraitPvalue_res.csv"))) 
  write.csv(moduleTraitCor_res, file = file.path(res_dir, paste0(tag_name, "_moduleTraitCor_res.csv")))
  
  # Create labelled heatmap
  
  data2return <- list(
    "moduleTraitCor_res" = moduleTraitCor_res,
    "moduleTraitPvalue_res" = moduleTraitPvalue_res,
    "MEs0_res"=MEs0_res
  )
  
  return(data2return)
}

wgcna_labeledHeatmap <- function(MEs0_res, moduleTraitCor_res, moduleTraitPvalue_res, sample_sheet_given, plot_dir, tag_name ) {
  library(WGCNA)
  
  #Print correlation heatmap between modules and traits
  
  #----------------------
  ## Labelled heatmap for all traits
  #----------------------
  textMatrix= paste(signif(moduleTraitCor_res, 2), "\n(",
                    signif(moduleTraitPvalue_res, 1), ")", sep= "")
  
  dim(textMatrix)= dim(moduleTraitCor_res)
  
  pdf(file.path(plot_dir, paste0(tag_name, "_labeledHeatmap.pdf")), paper = "A4r", width = 35, height = 12)
  
  labeledHeatmap(Matrix= moduleTraitCor_res,
                 xLabels= names(sample_sheet_given),
                 yLabels= names(MEs0_res),
                 ySymbols= names(MEs0_res),
                 colorLabels= FALSE,
                 colors= blueWhiteRed(50),
                 textMatrix= textMatrix,
                 setStdMargins= FALSE,
                 cex.text= 0.5,
                 zlim= c(-1,1),
                 main= paste0("Module-trait relationships:", tag_name, " mergeCutHeight:0.25"))
  dev.off()
  #----------------------
  
  #----------------------
  # Split
  #----------------------
  list_of_cols <- sample_sheet_splitter(sample_sheet_given)
  print(list_of_cols)
  
  for (i in names(list_of_cols)) {
    print(i)
    #Print correlation heatmap between modules and traits
    textMatrix= paste(signif(moduleTraitCor_res[,list_of_cols[[i]]], 2), "\n(",
                      signif(moduleTraitPvalue_res[,list_of_cols[[i]]], 1), ")", sep= "")
    
    dim(textMatrix)= dim(moduleTraitCor_res[,list_of_cols[[i]]])
    
    pdf(file.path(plot_dir, paste0(tag_name, "_", i, "_labeledHeatmap.pdf")), paper = "A4r", width = 35, height = 12)
    labeledHeatmap(Matrix= moduleTraitCor_res[,list_of_cols[[i]]],
                   xLabels= names(sample_sheet_given[,list_of_cols[[i]]]),
                   yLabels= names(MEs0_res),
                   ySymbols= names(MEs0_res),
                   colorLabels= FALSE,
                   colors= blueWhiteRed(50),
                   textMatrix= textMatrix,
                   setStdMargins= FALSE,
                   cex.text= 0.5,
                   zlim= c(-1,1),
                   main= paste0("Module-trait relationships:", tag_name, " mergeCutHeight:0.25"))
    dev.off()
  }
  
}

wgcna_dynamicModules <- function(picked_power, input_data.mat, tag_name, plot_dir) {
  library(WGCNA)
  
  #---------------------
  ## Dynamic tree cut: minClusterSize = 30
  #---------------------
  adjacency = adjacency(input_data.mat, power = picked_power)
  TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
  dissTOM = 1-TOM
  
  # Plot gene tree
  geneTree = hclust(as.dist(dissTOM), method = "average");
  
  pdf(file.path(plot_dir, paste0(tag_name, "_dynamic_tree.pdf")), paper = "A4r", width = 35, height = 12)
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);
  dev.off()
  
  # Module identification using dynamic tree cut
  # We like large modules, so we set the minimum module size relatively high:
  # minModuleSize = 30;
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                              pamRespectsDendro = FALSE, minClusterSize = 30);
  table(dynamicMods)
  length(table(dynamicMods)) 
  
  # Convert numeric labels into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  # Plot the dendrogram and colors underneath
  pdf(file.path(plot_dir, paste0(tag_name, "_dynamic_plotDendroAndColors.pdf")), paper = "A4r", width = 35, height = 12)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE, 
                      guideHang = 0.05,main = "Dynamic: Gene dendrogram and module colors")
  dev.off()
  #---------------------
  
  #---------------------
  ## compare with merge dynamic
  #---------------------
  MEList = moduleEigengenes(input_data.mat, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  
  # Merge close modules
  MEDissThres=0.40
  abline(h=MEDissThres, col = "red")
  merge = mergeCloseModules(input_data.mat, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
  mergedColors = merge$colors  
  mergedMEs = merge$newMEs  
  # Plot merged module tree
  pdf(file.path(plot_dir, paste0(tag_name, "_merge_dynamic_plotDendroAndColors.pdf")), paper = "A4r", width = 35, height = 12)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                      c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
  dev.off()
  #---------------------
  
  data2return <- list(
    "dynamicColors"=dynamicColors,
    "mergedColors"=mergedColors,
    "dissTOM"=dissTOM,
    "geneTree"=geneTree
  )
  
  return(data2return)
}
