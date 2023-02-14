#' @title
#' Get topologically associated domains (TADs)
#'
#' @description
#' Retrieve 3D genome browser TAD annotations. Downloaded on 22JUN08 from
#' http://3dgenome.fsm.northwestern.edu/downloads/hg38.TADs.zip.
#'
#' @param cell_types A vector of cell types to be returned. If not provided,
#' all types are returned.
#'
#' @import dplyr
#' @import magrittr
#' @import data.table
#' @import GenomicRanges
#'
#' @export

get_tads  <- function(cell_types){
  tad_file<- system.file("extdata","combined_tads_3dGB_2022.tsv",package="clugPac")
  #tad_file <- "/N/u/aclugston/resources/tads/hg38/combined_tads_3dGB_2022.tsv"
  tad_tb  <- fread(tad_file) %>%
    as_tibble %>%
    dplyr::mutate(cell_type = factor(cell_type)) %>%
    dplyr::filter(end > 0) #Some TAD annotations have a start coord but a zero end coord...not sure why but drop them.

  if(!missing(cell_types)){
    tad_tb<- dplyr::filter(tad_tb,cell_type %in% cell_types)
  }
  tad_tb %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    return()
}

