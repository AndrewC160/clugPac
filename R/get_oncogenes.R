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
#' "CancerGeneListSPCG_2023_formatted.tab" provided by Henry Martell. If set to
#' "oncokb"," this table is from the "OncoKB" database
#' (https://www.oncokb.org/) from the "oncokb_cancerGeneList.tsv" file.
#'
#' @param gene_source String of oncogene source table name. Defaults to "asclab".
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
  }else if(gene_source == "oncokb"){
    tb_out <- system.file("extdata","25JAN30-oncokb_cancer_gene_list.tsv",package="clugPac") %>%
      fread(header=TRUE) %>%
      as_tibble %>%
      dplyr::rename_all(function(x) gsub("[ -]","_",gsub("#","num",gsub("[()]","",tolower(x))))) %>%
      dplyr::rename(gene_name = hugo_symbol)

    bool_cols <- tb_out[1,] %>% unlist %>% sapply(function(x) x %in% c("Yes","No")) %>% which
    tb_out <- mutate_at(tb_out,bool_cols, function(x) ifelse(x == "Yes",TRUE,FALSE))
  }else{
    tb_out  <- system.file("extdata","CancerGeneListSPCG_2023_formatted.tab",package="clugPac") %>%
      fread(header=TRUE) %>%
      dplyr::rename(ensembl_id = ens,
                    gene_name = gene) %>%
      as_tibble
  }
  return(tb_out)
}
