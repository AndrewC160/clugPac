#' @title
#' Get oncogenes
#'
#' @description
#' Retrieve table of oncogenes and their IDs. If "source" is set to "bushman",
#' this table is retrieved from the Bushman Lab group's compilation of cancer
#' genes, downloaded from
#' http://www.bushmanlab.org/assets/doc/allOnco_June2021.tsv (Website:
#' http://www.bushmanlab.org/links/genelists). If set to
#' "asclab" (default), this table is from the file
#' "CancerGeneListSPCG_2023_formatted.tab" provided by Henry Martell.
#'
#'
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import data.table
#'
#' @export

get_oncogenes   <- function(gene_source="asclab"){
  if(gene_source == "bushman"){
    tb_out  <- system.file("extdata","all_oncogenes_June2021_formatted.tsv",package="clugPac") %>%
      fread(header=TRUE) %>%
      as_tibble
  }else{
    tb_out  <- system.file("extdata","CancerGeneListSPCG_2023_formatted.tab",package="clugPac") %>%
      fread(header=TRUE) %>%
      as_tibble
  }
  return(tb_out)
}
