#' @title 
#' Get exons.
#' 
#' @description 
#' Retrieve a table of exons annotated with Ensembl IDs. Exon table was 
#' generated using TxDB and biomaRt via the script:
#' 
#' '/N/u/aclugston/resources/gene_annotations/generate_exon_table.R.
#' 
#' @import data.table
#' @import tibble
#' @import magrittr
#' 
#' @export

get_exons  <- function(){
  exon_file <- system.file("extdata","exons.bed.gz",package="clugPac")
  
  fread(exon_file,sep="\t",header=FALSE,col.names=c("seqnames","start","end","width","strand","exon_id","ens_id")) %>%
    as_tibble %>%
    return()
}