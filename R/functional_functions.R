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
  
  CGP_curated_gene_sets <- subset(pathway_gene_sets, gs_subcat=="CGP")
  
  ## NABA, SA, SIG, WNT_
  CP_curated_gene_sets <- as.data.frame(subset(pathway_gene_sets, gs_subcat=="CP"))
  
  CP_Biocarta_curated_gene_sets <- as.data.frame(subset(pathway_gene_sets, gs_subcat=="CP:BIOCARTA"))
  CP_KEGG_curated_gene_sets <- as.data.frame(subset(pathway_gene_sets, gs_subcat=="CP:KEGG"))
  CP_PID_curated_gene_sets <- as.data.frame(subset(pathway_gene_sets, gs_subcat=="CP:PID"))
  CP_REACTOME_curated_gene_sets <- as.data.frame(subset(pathway_gene_sets, gs_subcat=="CP:REACTOME"))
  CP_WIKI_curated_gene_sets <- as.data.frame(subset(pathway_gene_sets, gs_subcat=="CP:WIKIPATHWAYS"))
  
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
  
  CGP_curated_gene_sets_list = split(x=CGP_curated_gene_sets$gene_symbol, 
                                     f=CGP_curated_gene_sets$gs_name)
  
  CP_curated_gene_sets_list = split(x=CP_curated_gene_sets$gene_symbol, 
                                    f=CP_curated_gene_sets$gs_name)
  
  CP_Biocarta_curated_gene_sets_list = split(x=CP_Biocarta_curated_gene_sets$gene_symbol, 
                                             f=CP_Biocarta_curated_gene_sets$gs_name)
  CP_KEGG_curated_gene_sets_list = split(x=CP_KEGG_curated_gene_sets$gene_symbol, 
                                         f=CP_KEGG_curated_gene_sets$gs_name)
  CP_REACTOME_curated_gene_sets_list = split(x=CP_REACTOME_curated_gene_sets$gene_symbol, 
                                             f=CP_REACTOME_curated_gene_sets$gs_name)
  CP_PID_curated_gene_sets_list = split(x=CP_PID_curated_gene_sets$gene_symbol, 
                                        f=CP_PID_curated_gene_sets$gs_name)
  CP_WIKI_curated_gene_sets_list = split(x=CP_WIKI_curated_gene_sets$gene_symbol, 
                                         f=CP_WIKI_curated_gene_sets$gs_name)
  
  
  
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
    "C2_CGP_curated" = CGP_curated_gene_sets_list,
    "C2_CP_curated_left"=CP_curated_gene_sets_list,
    "C2_BIOCARTA"=CP_Biocarta_curated_gene_sets_list,
    "C2_KEGG"=CP_KEGG_curated_gene_sets_list,
    "C2_REACTOME"=CP_REACTOME_curated_gene_sets_list,
    "C2_PID"=CP_PID_curated_gene_sets_list,
    "C2_WIKI"=CP_WIKI_curated_gene_sets_list,
    
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
rank_list_by <- function(table_data, option_given="logFC", GENE_SYMBOL.col="GENE_SYMBOL") {
  
  print(dim(table_data))
  table_data <-subset(table_data, !GENE_SYMBOL.col=="")
  print(table_data)
  print(dim(table_data))
  
  gene_list <- table_data[[option_given]]
  
  names(gene_list) = table_data[[GENE_SYMBOL.col]]
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
FGSEA_GSEA <- function(gene_list_provided, myGeneSet, title_given="example", nproc_given=2, minSize=10, maxSize=600, nPermSimple = 10000, padj.thres=0.25) {
  
  library(fgsea)
  fgRes <- fgsea::fgseaMultilevel(pathways = myGeneSet, minSize=minSize, maxSize=maxSize, nPermSimple = nPermSimple,
                                  stats = gene_list_provided, nproc = nproc_given) %>% as.data.frame()
  #minSize=1,
  #minSize=15, maxSize=600, nperm=10000) %>% as.data.frame()
  
  fgRes1 <- fgRes[fgRes$padj<padj.thres,]
  data.table::setorder(fgRes1, padj)

  returnList = list(
    fgRes = fgRes,
    p = list(
      "top" = plot_GSEA(head(fgRes1, n=10), title_given),
      "all" = plot_GSEA(fgRes1, title_given)
      )
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
  
  library(ggplot2)
  g = ggplot(fgRes, aes(reorder(pathway, ES), ES)) +
    geom_segment( aes(reorder(pathway, ES), xend=pathway, y=0, yend=ES)) +
    geom_point(aes( fill = Enrichment, size=padj),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                                 "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=title_given) + 
    theme_minimal()
  
  return(g)
}

#' Create GSEA loop
#' 
#' Plots GSEA results using ggplot and enrichment score provided
#' @param table_annot table to get results
#' @param folder_out folder to store results
#' @param name_given name for the comparison
#' @param nproc_given threads to use
#' @param dataSet.list list of GSEA datasets
#' @param GENE_SYMBOL.col name of the hgnc column
#' @param ranker.list columns to use for ranking
#' @param padj.thres Pvalue adjusted threshold
#' @export
FGSEA_GSEA_loop <- function(table_annot, folder_out, name_given, nproc_given, 
                            dataSet.list, GENE_SYMBOL.col="GENE_SYMBOL",  
                            ranker.list=c("log2FoldChange", "signed.logPval"), padj.thres=0.25, create.signed.logPval=FALSE) {
  
  library(openxlsx)
  
  dir.create(folder_out, showWarnings = FALSE)
  
  
  if (create.signed.logPval) {
    table_annot['signed.logPval'] <- sign(table_annot$log2FoldChange) * -log(table_annot$pvalue)  
  }
  
  ## FGSEA analysis
  fgsea2return = list()
  
  
  for (ranker in ranker.list ) {
    ## 
    print(paste0("+ Get data ranked by: ", ranker))
    gene_list_ranked <- rank_list_by(table_data = table_annot, 
                                     option_given = ranker, GENE_SYMBOL.col = GENE_SYMBOL.col)
    
    #print(gene_list_ranked)
    
    print("+ FGSEA analysis started...")
    print("")
    
    wb <- createWorkbook()
    xlsx_file <- file.path(folder_out, paste0(name_given, "-rank-by_", as.character(ranker), ".xlsx"))
    
    for (dataSet in names(dataSet.list)) {
      
      print(paste0("+ Analysis for: ", dataSet))
      FGSEA_data <- FGSEA_GSEA(na.omit(gene_list_ranked), dataSet.list[[dataSet]], 
                               title_given = paste0(name_given, 
                                                    " ranked by ", 
                                                    as.character(ranker), 
                                                    ":: Gene set = ", 
                                                    dataSet), 
                               nproc_given=nproc_given, padj.thres = padj.thres)
      
      file_name = paste0(name_given, "-rank-by_", as.character(ranker), "_", dataSet)
      # print data
      
      if (dim(subset(FGSEA_data$fgRes, padj<padj.thres))[1] > 1) {
        print("Printing plot in PDF..")
        save_pdf(folder_path = folder_out, name_file = paste0(file_name, ".top"), 
                 plot_given = FGSEA_data$p$top)
        
        save_pdf(folder_path = folder_out, 
                 name_file = paste0(file_name, ".all"), 
                 plot_given = FGSEA_data$p$all)
        
        addWorksheet(wb, dataSet)
        writeData(wb, dataSet, 
                  FGSEA_data$fgRes,
                  startRow = 2, startCol=2,rowNames = FALSE,
                  keepNA=TRUE,na.string="NA")
        
      }
      # 
      print("Printing data in XLSX..")
      print(xlsx_file)
      
      fgsea2return[[ranker]][[dataSet]] = FGSEA_data
      saveWorkbook(wb, xlsx_file, overwrite = TRUE)
    }
    
    
  }
  
  save(fgsea2return, file=file.path(folder_out, "gsea.RData"))
  return(fgsea2return)
  
}


#' Get gene annotation from BioMart
#'
#' Gets annotation for the set of genes desired
#' @param GeneList Gene IDs entries. ENSMBL genes only.
#' @param datasets By default: hsapiens_gene_ensembl
#' @export
get_gene_annotation <- function(Genelist, species="hsapiens_gene_ensembl", get_uniprot=FALSE) {
  ##get gene symbols from biomart - watch out for suffixes effect on annotation retrieval!!!
  
  library(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = species)
  
  ## If we filter by hgnc_symbol it gets duplicates and it is not uniq
  ## ENSEMBL id is unique
  
  if (get_uniprot) {
  
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
  } else {
    
    ## Just name, description and biotype    
    GeneSymbolsTable_full <-getBM(attributes = c('ensembl_gene_id', 
                                                 'hgnc_symbol', 'description', 'gene_biotype'),
                                  filters = 'ensembl_gene_id', 
                                  values = Genelist,
                                  mart = mart)
    
    
    return(GeneSymbolsTable_full)
  }
  
  
  
  
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


#' Create gene ontology enrichment for a given dataset
#'
#' Gets coordinates for the set of genes desired
#' @param list_of_sets List of dataset to use
#' @param geneUniverse List of total genes to use as "Universe"
#' @param GO_folder Folder to store results
#' @param tag2use Tag name to include
#' @param org.Db2use Database to use: Default: 'org.Hs.eg.db'. It should be loaded first: library('org.Hs.eg.db')
#' @param ont_types  List of ontologies classes to use. Default: c('CC', 'BP', 'MF') 
#' @param keyType.given SYMBOL as default. See keytypes(org.Hs.eg.db) as an example for putative keys to use
#' @export
enricher_loop <- function (list_of_sets, geneUniverse, GO_folder, tag2use,
                           org.Db2use='org.Hs.eg.db', ont_types = c('CC', 'BP', 'MF'), keyType.given='SYMBOL'){
  library(clusterProfiler)
  library(AnnotationDbi)
  library(openxlsx)
  
  wb <- createWorkbook()
  xlsx_file <- file.path(GO_folder, paste0(tag2use, "_enrichGO.xlsx"))
  
  
  results2save <-list()
  for (ont_T in ont_types) {
    for (set2test in names(list_of_sets)) {
      print(paste0("##### ", ont_T, ": ", set2test))
      
      name2use <- paste0(set2test, "-", ont_T)
      
      if (length(list_of_sets[[set2test]])<10) {
        print("enrichGO analysis: not possible. Too few genes provided")
        res <- ""
        p <- ""
        dat <- ""
        
      } else {
        print("enrichGO analysis")
        res <- clusterProfiler::enrichGO(gene = list_of_sets[[set2test]], 
                                         universe = geneUniverse,
                                         OrgDb = org.Db2use, 
                                         keyType = keyType.given, 
                                         ont = ont_T)
        
        
        print("save results")
        dat <- head(res, n = 100)
        
        if (dim(dat)[1]>10) {
          p <- enrichplot::dotplot(res)
          save_pdf(p, folder_path = GO_folder, name_file = paste0(tag2use, "_", name2use)  )
        }
        # save results
        write.csv(dat, file=file.path(GO_folder, paste0(tag2use, "_", name2use, ".csv")))
        
        addWorksheet(wb, name2use)
        writeData(wb, name2use, dat,
                  startRow = 2, startCol=2,rowNames = FALSE,
                  keepNA=TRUE,na.string="NA")
        
      }
      
      #
      results2save[[set2test]][[ont_T]] <- list(
        'results' = res,
        'dat' = dat,
        'plot' = p
      )
    }
  }
  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  
  return(results2save)  
}

