#' Create limma DE analysis 
#'
#' This functions creates limma DE analysis
#' @param normData dataframe to use as normalized values 
#' @param design.given Matrix with design 
#' @param contrasts.matrix.given makeContrast object to test
#' @param coef.given List of coefficients names to test.
#' @export

limma_DE_function <- function(normData, design.given, contrasts.matrix.given, coef.given) {
  # model.matrix: 
  # makeContrasts: make all-pairwise comparisons between the gruoups
  
  # lmFit: estimate the fold changes and standard errors by fitting a linear model for each gene.
  # eBayes: apply empirical Bayes smoothing to the standard errors
  # topTable: show statistical results
  # contrasts.fit: make all-pairwise comparisons between the gruoups
  
  fit_1 <- limma::lmFit(normData, design.given)
  fit_2 <- limma::contrasts.fit(fit_1, contrasts.matrix.given)
  efit_3 <- limma::eBayes(fit_2, trend=TRUE)        # empirical Bayes adjustment
  results_fit <- decideTests(efit_3)
  
  results_list <- list()
  for (c in coef.given) {
    all_data.res = limma::topTable(efit_3, number = 3e6, 
                                   adjust.method = "BH", coef = c)
    print(head(all_data.res))
    results_list[[c]][['all.data']] = all_data.res
    
    filt.res <- filter_signficant_limma(all_data.res)
    results_list[[c]][['sign.data']] = filt.res
    results_list[[c]][['sign.genes']] = rownames(filt.res)
    results_list[[c]][['sign.count']] = length(filt.res)
    
  }
  
  results2return <- list(
    "vennDiagram" = vennDiagram(results_fit),
    "results_list" = results_list
    
  )
  
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



