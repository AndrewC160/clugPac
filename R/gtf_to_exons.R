#' @title
#' GTF to exons
#'
#' @description
#' Given a GTF file, extract all exons within region(s) specified by <granges_query>.
#' Returns a GRanges object. If no GTF file is provided, data is retrieved from
#' packaged "Homo_sapiens.GRCh38.104.chr.tabix.gtf.gz" file.
#'
#' @param granges_query GenomicRanges object of genomic regions to be searched.
#' @param tabix_gtf Tabix-indexed GTF file.
#' @param cache_tsv TSV file in which to cache exon results, if not provided no cache will be written/read.
#' @param overwrite_cache Should existing cache files be overwritten? Defaults to FALSE.
#' @param as_grange Should data be returned as a GRanges object? Defaults to TRUE, otherwise is returned as a tibble.
#'
#' @import data.table
#' @import dplyr
#' @import tidyr
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomeInfoDb
#' @import stringr
#' @import magrittr
#'
#' @export

gtf_to_exons    <- function(granges_query,tabix_gtf,cache_tsv=NULL,overwrite_cache=FALSE,as_grange=TRUE){
  select  <- dplyr::select
  filter  <- dplyr::filter
  mutate  <- dplyr::mutate

  if(missing(tabix_gtf)){
    tabix_gtf   <- system.file("extdata","Homo_sapiens.GRCh38.104.chr.tabix.gtf.gz",package="clugPac")
  }
  if(!is.null(cache_tsv) & !overwrite_cache){
    if(file.exists(cache_tsv)){
      tb_exons  <- fread(cache_tsv,sep="\t",header = TRUE)
    }else{
      tb_exons  <- NULL
    }
  }else{
    tb_exons  <- NULL
  }
  if(is.null(tb_exons)){
    seqlevelsStyle(granges_query) <- "NCBI"
    txt <- Rsamtools::scanTabix(tabix_gtf,param = granges_query)
    tb_exons  <- lapply(names(txt), function(text_table){
      fread(text=txt[[text_table]],
            sep="\t",
            header=FALSE,
            col.names = c("seqnames",
                          "source",
                          "type",
                          "start",
                          "end",
                          "V6",
                          "strand",
                          "V8",
                          "vals")) %>%
        as_tibble %>%
        select(-V6,-V8) %>%
        filter(type == "exon") %>%
        mutate(seqnames= paste0("chr",seqnames),
               query_range=text_table)
    }) %>%
      do.call(rbind,.) %>%
      as_tibble %>%
      rowwise %>%
      mutate(vals = list(unlist(str_split(vals,pattern = ";")))) %>%
      ungroup %>%
      unnest(vals) %>%
      mutate(vals = trimws(vals),
             vals = gsub('"',"",vals),
             name = str_extract(vals,"^[[:alpha:]_]+"),
             val = str_extract(vals,"[^ ]+$")) %>%
      select(-vals) %>%
      filter(!is.na(name)) %>%
      filter(name %in% c("gene_name",
                         "gene_id",
                         "transcript_id",
                         "ens_id",
                         "exon_number",
                         "exon_id")) %>%
      distinct %>%
      pivot_wider(id_cols = c(seqnames,start,end,strand,query_range),
                  names_from = name,values_from = val,
                  values_fn = function(x) list(unique(x)))
    #Some exons don't have gene names, so in the off chance this happens in ALL exons one needs to add the gene_name column manually.
    if(!"gene_name" %in% tb_exons) tb_exons <- mutate(tb_exons,gene_name = list(""))
    tb_exons <- tb_exons %>%
      unnest(exon_number) %>%
      unnest(exon_id) %>%
      unnest(gene_name) %>%
      unnest(gene_id) %>%
      unnest(transcript_id) %>%
      select(gene_name,gene_id,transcript_id,exon_id,exon_number,everything(),query_range)
    if(!is.null(cache_tsv)){
      write.table(tb_exons,file = cache_tsv,quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")
    }
  }
  if(as_grange){
    tb_exons <- makeGRangesFromDataFrame(tb_exons,keep.extra.columns = TRUE)
  }
  return(tb_exons)
}
