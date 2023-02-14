#' @title
#' Get oncogenes
#'
#' @description
#' Retrieve the Bushman Lab group's compilation of cancer genes, downloaded
#' from http://www.bushmanlab.org/assets/doc/allOnco_June2021.tsv.
#'
#' Website: http://www.bushmanlab.org/links/genelists
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import data.table
#'
#' @export

get_oncogenes   <- function(){
  gene_file <- system.file("extdata","all_oncogenes_June2021_formatted.tsv",package="clugPac")
  fread(gene_file,header=TRUE) %>%
    as_tibble %>%
    return()
}
