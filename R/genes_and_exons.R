#' @title
#' Get genes and exon data.
#'
#' @description
#' Retrive gene and exon data from a GTF file and format as GRanges if
#' <as_granges> is TRUE. If no Grange region is provided, all genes are
#' retrieved, but exons are not.
#'
#' @param grange_win GRange denoting region within which to find genes and exons.
#' @param return_exons Should exons be retrieved as well? Defaults to FALSE.
#' @param clip_features Should features be truncated to fit strictly within the window specified? Defaults to FALSE.
#' @param cache_gene_tsv TSV file in which to cache genes via clugPac::gtf_to_genes(). If not provided, no cache is created/read.
#' @param cache_exon_tsv TSV file in which to cache exons via clugPac::gtf_to_exons(). If not provided, no cache is created/read.
#' @param overwrite_cache Should existing cache TSVs be overwritten if found? Defaults to FALSE.
#'
#' @import GenomicRanges
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @import stringr
#'
#' @export

genes_and_exons <- function(grange_win,return_exons=TRUE,clip_features = TRUE,cache_gene_tsv=NULL,cache_exon_tsv=NULL,overwrite_cache = FALSE,as_grange = FALSE){
  tb_genes    <- gtf_to_genes(cache_file = cache_gene_tsv,overwrite_cache = overwrite_cache,as_grange = FALSE)
  tb_exons    <- NULL
  onco_ens_ids<- get_oncogenes()$ensembl_id

  #If a window of interest is provided, clip genes to that window.
  #Also, get exon data. Do Not get exon data otherwise (too many to be reasonable).
  if(!missing(grange_win)){
    if(clip_features){
      tb_genes<- clip_granges(grange_in = tb_genes,grange_window = grange_win,replace_cols = TRUE,include_clipped_col = TRUE)
    }
    if(return_exons){
      #Get exons, and get as GRanges if needed.
      tb_exons<- gtf_to_exons(granges_query = grange_win,cache_tsv = cache_exon_tsv,overwrite_cache = overwrite_cache,as_grange = as_grange)
      if(clip_features){
        tb_exons  <- clip_granges(tb_exons,grange_window = grange_win,replace_cols = TRUE,include_clipped_col = TRUE)
      }
      #Label oncogene exons.
      tb_exons$bushman_onco <- tb_exons$gene_id %in% onco_ens_ids
    }
  }else if(return_exons){
    warning("genes_and_exons() does not return exons for whole genome, use gtf_to_exons() for large queries.")
  }
  #Label oncogenes.
  tb_genes<- mutate(tb_genes,bushman_onco = gene_id %in% onco_ens_ids)

  #Convert gene tibble to GRanges, if needed.
  if(as_grange){
    tb_genes  <- makeGRangesFromDataFrame(tb_genes,keep.extra.columns = TRUE)
  }
  if(is.null(tb_exons)){
    return(tb_genes)
  }else{
    return(list(genes=tb_genes,exons=tb_exons))
  }
}
