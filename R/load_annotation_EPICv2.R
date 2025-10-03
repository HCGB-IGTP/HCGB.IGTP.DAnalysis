#' Load different objects into a GenomicRanges list
#'
#' @param EPIC_annotv2 
#' @param re_annotated 
#' @param elements_UCSC 
#' @param gtf_file 
#'
#' @export
 load_annotation_EPICv2 <- function(
    EPIC_annotv2 = "/imppc/labs/lslab/share/data/references/EPIC/EPIC_v2/EPICv2.annot.RData",
    re_annotated = "/imppc/labs/lslab/share/data/references/EPIC/EPIC_v2/EPICv2_reannotated_manifest_v2.0.csv",
    elements_UCSC = "/imppc/labs/lslab/share/data/references/EPIC/EPIC_v2/elements_UCSC_annot.RData",
    gtf_file = "/imppc/labs/lslab/share/data/references/Human_reference_latest/Homo_sapiens.GRCh38.108.chr._ncRNA-RNAcentral.gtf") {

  ###################################
  ## get annotation
  ###################################
  
  load(EPIC_annotv2)
  
  EPIC_MEsteller$Probe_name  <- stringr::str_remove(EPIC_MEsteller$ProbeID, "_.*")
  
  EPIC_MEsteller.GR <- GenomicRanges::makeGRangesFromDataFrame(EPIC_MEsteller,  
                                                               start.field = "CpGbeg", end.field = "CpGend", 
                                                               seqnames.field =  "CpGchrm", 
                                                               strand.field = "probeStrand",
                                                               keep.extra.columns = TRUE)
  
  ## Later if necessary load other annotations
  re_annotated.df <- read.csv(re_annotated)
  head(re_annotated.df)
  rownames(re_annotated.df) <- re_annotated.df$IlmnID
  
  re_annotated.GR <- GenomicRanges::makeGRangesFromDataFrame(re_annotated.df,  
                                                             start.field = "MAPINFO", end.field = "MAPINFO", 
                                                             seqnames.field =  "CHR",
                                                             keep.extra.columns = TRUE, na.rm = TRUE)
  
  ## load annotation
  load(file = elements_UCSC)
  # RefseqFuncElemns_gtf
  # encodeCcreCombined_gtf
  
  ## genome annotation
  gtf_GRCh38 <- rtracklayer::import(gtf_file)
  
  list_annot.GR <- list(
    "EPIC_MEsteller.GR" = EPIC_MEsteller.GR, 
    "re_annotated.GR" = re_annotated.GR,
    "RefseqFuncElemns_gtf" = RefseqFuncElemns_gtf,
    "encodeCcreCombined_gtf" = encodeCcreCombined_gtf,
    "GTF_genome" = gtf_GRCh38)
  
  ##-------------------------------
  ## Fix seqnames not well annotated
  ##-------------------------------
  seqnames(EPIC_MEsteller.GR)
  seqlevelsStyle(list_annot.GR$GTF_genome) <- 'NCBI'
  
  gr1 = list_annot.GR$GTF_genome[seqnames(list_annot.GR$GTF_genome) %in% c(1:22, "X", "Y")]
  seqlevels(gr1) = as.character(unique(seqnames(gr1)))
  seqlevels(gr1)
  seqlevels(gr1) <- paste0("chr", seqlevels(gr1))
  
  list_annot.GR$GTF_genome <- gr1
  ##-------------------------------
  
  list_annot.GR
  
  ###############################
}
