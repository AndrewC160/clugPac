#' @title Transcription start site distances
#'
#' @description
#' Measure distances between a set of genomic regions and their nearest TSS. Negative
#' values indicate upstream respective of gene strand. Data is returned as a vector of integers. A
#' column name of gene IDs can be provided via <gene_id_col> which will be used to name the output
#' vector.
#'
#' @param gr_in GenomicRanges object of regions to measure.
#' @param gr_genes_in GenomicRanges object of gene annotations to measure against.
#' @param gene_id_col Name of metacolumn for <gr_genes_in> which will be used as gene IDs. Values will be returned as names of distance vector.
#'
#' @import tidyr
#' @import dplyr
#' @import magrittr
#' @import GenomicRanges
#' @import tibble
#'
#' @export

tss_distances <- function(gr_in,gr_genes_in,gene_id_col){
  if(missing(gr_genes_in)) stop("GRanges of gene annotations must be provided.")
  gr_tss  <- GenomicRanges::promoters(gr_genes_in,upstream = 0,downstream = 0)

  gene_col  <- paste0("tss_",c(1:length(gr_tss)))
  if(!missing(gene_id_col)){
    if(!gene_id_col %in% colnames(mcols(gr_genes_in))){
      message("'",gene_id_col,"' not found in gr_genes_in, skipping IDs.")
    }else{
      gene_col <- mcols(gr_genes_in)[[gene_id_col]]
    }
  }
  mcols(gr_tss)$id <- gene_col

  o_laps  <- GenomicRanges::distanceToNearest(gr_in,gr_tss)

  as_tibble(o_laps) %>%
    mutate(tss_strand = as.character(GenomicRanges::strand(gr_tss))[subjectHits],
           tss_pos = GenomicRanges::start(gr_tss)[subjectHits],
           pk_start = GenomicRanges::start(gr_in)[queryHits],
           pk_end = GenomicRanges::end(gr_in)[queryHits],
           id = mcols(gr_tss)$id[queryHits]) %>%
    mutate(is_upstream = case_when(tss_strand == "+" & pk_end < tss_pos ~ TRUE,
                                   tss_strand == "-" & pk_start > tss_pos ~ TRUE,
                                   TRUE ~ FALSE)) %>%
    mutate(distance = ifelse(is_upstream,-distance,distance)) %>%
    vectify(distance,id) %>%
    return()
}
