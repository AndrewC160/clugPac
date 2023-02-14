#' @title
#' Get Ensembl IDs
#'
#' @description
#' Retrieve Ensembl IDs from specific genes, if provided, otherwise return a
#' full table. IDs are retrieved from the same GTF files used by gtf_to_genes.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import data.table
#'
#' @export

get_ensembl_ids <- function(gene_symbols){
  # if(is.null(id_cache_table)){
  #   id_cache_table  <- system.file("extdata","gene_name_ensembl_key.tsv",package="clugPac")
  # }
  # if(file.exists(id_cache_table)){
  #   tb_ids  <- fread(id_cache_table,sep="\t",header=TRUE) %>% as_tibble
  # }else{
  #   tb_ids   <- gtf_to_genes() %>%
  #     as_tibble %>%
  #     select(gene_name,gene_id) %>%
  #     arrange(gene_name) %T>%
  #       write.table(file = id_cache_table,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  # }
  id_cache_table  <- system.file("extdata","gene_name_ensembl_key.tsv",package="clugPac")
  tb_ids  <- fread(id_cache_table,sep="\t",header=TRUE) %>% as_tibble
  if(!missing(gene_symbols)){
    ids <- tb_ids %>%
      dplyr::filter(!gene_name=="") %>%
      vectify(gene_id,gene_name)
    out_p <- ids[gene_symbols] %>%
      setNames(gene_symbols)
  }else{
    out_p <- tb_ids
  }
  return(out_p)
}
