#' Create limma DE analysis for all coefficients 
#'
#' This functions creates limma DE analysis
#' @param normData dataframe to use as normalized values 
#' @param df_treatment_Ind Sample sheet with information for plotting
#' @param design.given Matrix with design 
#' @param contrasts.matrix.given makeContrast object to test
#' @param coef.given.list List of coefficients names to test.
#' @param OUTPUT_Data_sample Output directory to save results
#' @param comp_name Column name in sample sheet given to test.
#' @param pCutoff.given Adjusted pvlaue cutoff. Default=0.05, 
#' @param FCcutoff.given Log Fold change cutoff. Default=log2(1.2), 
#' @param forceResults Set flag to force results if previously generated. Default=FALSE
#' @param EPIC Boolean to indicate EPIC probes, avoid generating huge xlsx file. Default=FALSE
#' @export
get_limma_results <- function(normData, df_treatment_Ind, design.given, contrasts.matrix.given, 
                              coef.given.list, OUTPUT_Data_sample,
                              comp_name,  pCutoff.given=0.05, 
                              FCcutoff.given=log2(1.2), forceResults=FALSE, EPIC=FALSE, paired_variable=NULL) {
  
  
  # model.matrix: 
  # makeContrasts: make all-pairwise comparisons between the gruoups
  
  # lmFit: estimate the fold changes and standard errors by fitting a linear model for each gene.
  # eBayes: apply empirical Bayes smoothing to the standard errors
  # topTable: show statistical results
  # contrasts.fit: make all-pairwise comparisons between the gruoups
  
  corfit <- ""
  if (is.null(paired_variable)) {
    fit_1 <- limma::lmFit(normData, design.given)
    
  } else {
    
    ## estimate correlation between variables on the same subject
    corfit <- duplicateCorrelation(normData, design.given, block=paired_variable)
    fit_1 <- limma::lmFit(normData, design.given, block = paired_variable, correlation = corfit$consensus)  
  }
  
  fit_2 <- limma::contrasts.fit(fit_1, contrasts.matrix.given)
  efit_3 <- limma::eBayes(fit_2, trend=TRUE)        # empirical Bayes adjustment
  results_fit <- decideTests(efit_3)
  
  results_list <- list()
  for (c in coef.given.list) {
    
    print(paste0("Get results for coefficient: ", c))
    coeff_split = unlist(strsplit(c, "_vs_")[1])
    numerator = coeff_split[1]
    denominator = coeff_split[2]
    
    c.dir <- file.path(OUTPUT_Data_sample, c)
    dir.create(c.dir)
    
    results.limma <- limma_DE_function(efit_3 = efit_3, normData =  normData, 
                                       df_treatment_Ind = df_treatment_Ind, design.given = design.given, 
                                       contrasts.matrix.given = contrasts.matrix.given, 
                                       coef.given = c, OUTPUT_Data_sample = c.dir,
                                       comp_name = comp_name, 
                                       numerator=numerator, 
                                       denominator = denominator,
                                       pCutoff.given=pCutoff.given, 
                                       FCcutoff.given=FCcutoff.given, 
                                       forceResults=forceResults, EPIC=EPIC)
    
    results_list[[c]] = results.limma
    
  }
  
  results2return <- list(
    "efit"=efit_3,
    "results_fit" = results_fit,
    #"vennDiagram" = vennDiagram(results_fit),
    "results_list" = results_list,
    "normData" = normData,
    "design" = design.given,
    "contrasts.matrix" = contrasts.matrix.given,
    "coef.given" = coef.given.list
  )
  
  if (!is.null(paired_variable)) {
    results2return[['corrfit']] = corfit 
  }
  
  return(results2return)
}

#' Create limma DE analysis 
#'
#' This functions creates limma DE analysis
#' @param efit_3 Ebayes fit results
#' @param normData dataframe to use as normalized values 
#' @param df_treatment_Ind Sample sheet with information for plotting
#' @param design.given Matrix with design 
#' @param contrasts.matrix.given makeContrast object to test
#' @param coef.given One coefficient to test.
#' @param OUTPUT_Data_sample Output directory to save results
#' @param comp_name Column name in sample sheet given to test.
#' @param numerator Comparison numerator.
#' @param denominator Comparison denominator.
#' @param pCutoff.given Adjusted pvlaue cutoff. Default=0.05, 
#' @param FCcutoff.given Log Fold change cutoff. Default=log2(1.2), 
#' @param forceResults Set flag to force results if previously generated. Default=FALSE
#' @param EPIC Boolean to indicate EPIC probes, avoid generating huge xlsx file. Default=FALSE

#' @export
limma_DE_function <- function(efit_3, normData, df_treatment_Ind, design.given, contrasts.matrix.given, 
                              coef.given, OUTPUT_Data_sample,
                              comp_name, numerator, denominator,
                              pCutoff.given=0.05, FCcutoff.given=log2(1.2), forceResults=FALSE, EPIC=FALSE) {
  
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(pheatmap)
  library(openxlsx)
  
  if (forceResults) {
    
  } else {
    ## check if previously done
    if (file.exists(file.path(OUTPUT_Data_sample, "data2return.RData"))) {
      print("Data already available in:")
      print(OUTPUT_Data_sample)
      return()
    }
  }
  
  #--------------------------
  ## Prepare data
  #--------------------------
  results_list <- list()
  all_data.res = limma::topTable(efit_3, number = 3e6, 
                                 adjust.method = "BH", coef = coef.given)
  print(head(all_data.res))
  
  ## Merge normalized values and differential expression
  alldata.norm.res <- merge(normData, all_data.res, by="row.names", sort=FALSE)
  names(alldata.norm.res)[1] <- "Gene"
  
  # filter
  filt.res <- filter_signficant_limma(alldata.norm.res, sign_value = pCutoff.given, LFC=FCcutoff.given)
  names(filt.res)[1] <- "Gene"
  rownames(filt.res) <- filt.res$Gene
  filt.res$Gene <- NULL
  #--------------------------
  
  #--------------------------
  # Get only for this samples
  #--------------------------
  
  listOfSampls <- c()
  
  if (denominator=="reference") {
    
    listOfSampls <- c("None")
    file_name <- paste0("cmp.", coef.given)
    name.cmp <- comp_name
    Samplslist <- c("All")
    
  } else {
    print("Comparison: ")
    print(comp_name)
    
    # c: coefficient contains: XXX_vs_YYY
    
    print("numerator: ")
    print(numerator)
    
    print("denominator: ")
    print(denominator)
    
    print("Samples: ")
    subsheet <- df_treatment_Ind[comp_name]
    colnames(subsheet)[1] <- 'comp_name'
    print(rownames(subsheet))
    
    listOfSampls <- c(rownames(subset(subsheet, comp_name==numerator)),
                      rownames(subset(subsheet, comp_name==denominator)))
    
    Samplslist <- list(
      "comp_name" = comp_name,
      "numerator" = rownames(subset(subsheet, comp_name==numerator)),
      "denominator" = rownames(subset(subsheet, comp_name==denominator))
    )
    
    file_name <- paste0(comp_name, "_",  numerator, "_vs_", denominator)
    name.cmp <- paste0(comp_name, ": ",  numerator, " vs. ", denominator)
  }
  #--------------------------
  
  #--------------------------
  # Write Results tables
  #--------------------------
  if (EPIC) {
    print("No saving csv table with all results, just save as RData")
  } else {
    write.table(alldata.norm.res, 
                file=file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_NormValues-and-DE_table.txt")), 
                sep="\t", row.names=T, col.names=NA, quote=T)
  }  
  write.table(filt.res, file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting_table_SignificantDE.txt")), 
              sep="\t", row.names=T, col.names=NA, quote=TRUE)
  
  
  if (EPIC) {
    print("No XLSX file generated for EPIC results")
  } else {
    
    DE.filename <- file.path(OUTPUT_Data_sample, paste0(file_name, "-ResultsCounting.xlsx"))
    sheet_name <- file_name
    title_name <- paste0("Comparison for coefficient: ", coef.given)
    
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, title_name, startRow = 2, startCol=2,
                        rowNames = FALSE, keepNA=TRUE,na.string="NA")
    openxlsx::writeData(wb, sheet_name, alldata.norm.res, startRow = 4, startCol=2, 
                        rowNames = FALSE, keepNA=TRUE,na.string="NA")
    openxlsx::saveWorkbook(wb, DE.filename, overwrite = TRUE)
    
  }
  #--------------------------
    
  #--------------------------
  ## add volcano plot
  #--------------------------
  plt_volcano <- EnhancedVolcano::EnhancedVolcano(all_data.res, 
                                                  x="logFC", y="adj.P.Val", lab="",
                                                  pCutoff=pCutoff.given, FCcutoff=FCcutoff.given, pointSize=3, labSize=6) + 
    ggplot2::scale_x_continuous() + 
    ggplot2::labs(title = paste0("Comparison for ", name.cmp) )
  
  if (EPIC) {
    #print("No Volcano file generated for EPIC results")
    HCGB.IGTP.DAnalysis::save_png(folder_path = OUTPUT_Data_sample, 
                                name_file = paste0(file_name, "_DiffExpression-volcano-plot"), 
                                plot_given = plt_volcano)

  } else {
    HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                name_file = paste0(file_name, "_DiffExpression-volcano-plot"), 
                                plot_given = plt_volcano)

  }
  #--------------------------
  
  #--------------------------
  # Heatmap
  #--------------------------
  if (length(rownames(filt.res))>10) {
    
    print(head(filt.res))
    
    ## plot rld
    data2heatmap = filt.res[,!colnames(filt.res) %in% c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val","B")] %>% na.omit()
    print(data2heatmap[1:50, 1:10])
    #print(tail(colnames(data2heatmap)))
    
    plot1 <- ""    
    plot1 <- try(pheatmap(head(data2heatmap, 50), main="NormCounts Pheatmap (p.adj<0.05 and [FC]>1.2)",
                      cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                      annotation_col = df_treatment_Ind,
                      color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
                      scale="row" ## centered and scale values per row
    ))

    HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                  name_file = paste0(file_name, "_top50_DEgenes_Heatmap_allSamples"), plot_given = plot1)
    results_list[['heatmap.all']] = plot1
    
    
    if (listOfSampls[1]=="None") {
      print("")
      
    } else {
      ## Only samples included in comparison
      ## plot rld
      data2heatmap2 = data2heatmap[,c(Samplslist$numerator, Samplslist$denominator)]
      plot2 <- ""
      plot2 <- try(pheatmap(head(data2heatmap2, 50), main="NormCounts Pheatmap (p.adj<0.05 and [FC]>1.2)",
                        cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames = TRUE, legend = TRUE,
                        annotation_col = df_treatment_Ind,
                        color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)), 
                        scale="row" ## centered and scale values per row
      ))
      HCGB.IGTP.DAnalysis::save_pdf(folder_path = OUTPUT_Data_sample, 
                                    name_file = paste0(file_name, "_top50_DEgenes_Heatmap_Samples"), 
                                    plot_given = plot2)
      results_list[['heatmap.only']] = plot2
      
    }
    
  }
  #--------------------------
  
  #--------------------------
  # save results in list
  #--------------------------
  results_list[['all.data']] = all_data.res
  results_list[['alldata.norm.res']] = alldata.norm.res
  results_list[['sign.data']] = filt.res
  results_list[['sign.genes']] = rownames(filt.res)
  results_list[['sign.count']] = length(rownames(filt.res))
  results_list[['plt_volcano']] = plt_volcano
  results_list[['Samplslist']] = Samplslist
  #--------------------------
  
  ######################################################################
  print ("Finish here for: ")
  print(name.cmp)
  ######################################################################
  
  data2save <- list(
    "Samplslist" = Samplslist,
    "alldata" = all_data.res,
    "volcan_plot" = plt_volcano,
    "alldata.norm.res"=alldata.norm.res,
    "sign.df"=filt.res,
    "sign.genes"=rownames(filt.res),
    "sign.count"=length(rownames(filt.res))
  )
  
  ## dump in disk RData
  save(data2save, file=file.path(OUTPUT_Data_sample, "data2return.RData"))
  
  return(results_list)
}

#' Filter limma DE analysis 
#'
#' This functions filters limma DE results
#' @param normData results limma dataframe to use as normalized values 
#' @param sign_value Pvalue adjusted cutoff: Default=0.05
#' @param LFC Log Fold Change cutoff: Default: 0.26
#' @export
filter_signficant_limma <- function(dataF, sign_value = 0.05, LFC=0.26) {
  ## limma
  dataFilt <- subset(dataF, abs(logFC)>LFC & adj.P.Val<sign_value)
  dataFilt <- dataFilt[order(dataFilt$adj.P.Val),]
  return(dataFilt)
}


#' Make all pairwise contrast for limma
#' 
#' Original: https://bioinformatics.stackexchange.com/questions/18570/generating-contrast-matrix-for-limma-in-loop
#'
#' @param group Vector of names to use 
#' @param delim Character to use for delimiting: "_vs_" as default
#'
#' @export
make_all_contrasts <- function (group, delim="_vs_"){
  
  suppressMessages(require(limma))
  
  #/ ensure that group levels are unique
  group <- sort(unique(as.character(group)))
  
  #/ make all combinations
  cb   <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
  
  #/ make contrasts
  contrasts<- limma::makeContrasts(contrasts=cb, levels=group)
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  
  return(contrasts)
}


