## Get MEAL Data
library(MEAL)
library(MultiDataSet)
library(missMethyl)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(rexposome)

prepare_MEAL_data <- function(betas, counts, sample_sheet, list_samples, comp1, comp2) {
  
  ## subset comparison1
  list1 <- subset(sample_sheet, Study==comp1)$Sample_Name
  list1 <- list1[list1 %in% list_samples]
  
  print("+ Create multidataset for Comparison:")
  print(comp1)
  print("+ Use samples:")
  print(list1)
  
  multi_set1 <- get_MEAL_data(betas[,list1], counts[,list1], sample_sheet[list1,], list1)
  
  ## subset comparison2
  list2 <- subset(sample_sheet, Study==comp2)$Sample_Name
  list2 <- list2[list2 %in% list_samples]
  
  multi_set2 <- get_MEAL_data(betas[,list2], counts[,list2], sample_sheet[list2,], list2)
  
  list2return <- list(
    "comp1" = multi_set1,
    "comp2" = multi_set2
  )
  
  return(list2return)
  
  
}

###
get_MEAL_data <- function(betas, counts, sample_sheet, list_samples) {
  ## use annotation according to EPIC and human
  GRset <- minfi::makeGenomicRatioSetFromMatrix(mat=as.matrix(betas), 
                                                array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
  
  class(GRset)
  
  ## change sequence level style
  seqlevelsStyle(GRset)
  #[1] "UCSC"
  seqlevelsStyle(GRset) <- 'NCBI'
  
  
  ## get annotation
  #rownames(counts) <- counts$Gene
  coordinates_full <- get_gene_coordinates(rownames(counts))
  rownames(coordinates_full) <- coordinates_full$ensembl_gene_id
  colnames(coordinates_full)[5] <- 'old_strand'
  colnames(coordinates_full)[2] <- 'chromosome'
  colnames(coordinates_full)[3] <- 'start'
  colnames(coordinates_full)[4] <- 'end'
  
  coordinates_full$strand <- ifelse(coordinates_full$old_strand == '-1', '-', '+')
  head(coordinates_full)
  coordinates_full$old_strand <- NULL
  coordinates_full$ensembl_gene_id <- NULL
  
  coordinates_full <- coordinates_full[order(row.names(coordinates_full)),]
  
  print(head(coordinates_full))
  
  ## get annotated genes
  exprs_fixed <- counts[ rownames(counts) %in% rownames(coordinates_full),]
  dim(exprs_fixed)
  exprs_fixed <- exprs_fixed[, list_samples]
  dim(exprs_fixed)
  
  exprs_fixed <- exprs_fixed[order(row.names(exprs_fixed)),]
  exprs_fixed$Gene <- NULL
  exprs_fixed$baseMean <- NULL
  exprs_fixed$pvalue <- NULL
  exprs_fixed$padj <- NULL
  exprs_fixed$log2FoldChange <- NULL
  exprs_fixed$lfcSE <- NULL
  
  print(dim(exprs_fixed))
  
  ## get samples to compare
  phenoData1 <- new("AnnotatedDataFrame", data=sample_sheet)
  
  ## create expression set
  print (colnames(exprs_fixed))
  print (colnames(coordinates_full))
  print (rownames(phenoData1))
  
  expr_set1 <- ExpressionSet(assayData = as.matrix(exprs_fixed), 
                             phenoData = phenoData1)
  print (class(expr_set1))
  print (dim(expr_set1))
  print (head(expr_set1))
  
  ## add features
  fData(expr_set1) <- coordinates_full
  
  multi <- createMultiDataSet()
  multi <- add_genexp(multi, expr_set1)
  multi <- add_methy(multi, GRset)
  
  return(multi)
}


## analysis
multidataset <- prepare_MEAL_data(betas_df, exprs_df, targetsFile, samples_both, "comp1", "comp2")

## Filter using genomic range
targetRange <- GRanges('4:151500192-151510000')
test1 <- multidataset$multi_comp1[, , targetRange]
methExprs <- correlationMethExprs(test1)

## Filter using list of CPGs to test
list_cpgs <- c('', '')
correlationMethExprs(multiset = multidataset$comp1, 
                     sel_cpgs = list_cpgs)


