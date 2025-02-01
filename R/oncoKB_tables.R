#' OncoKB tables
#'
#' @description
#' Retrieve OncoKB data table(s) as retrieved on 25JAN30. Tables include gene-
#' specific therapies (therapies), recognized alterations (variants), oncogenes
#' (oncogenes), drug/variant level descriptions (levels), and FDA-approved drug
#' descriptions ("drug_associations").
#'
#' @param table_names Vector of table names to retrieve.
#'
#' @import dplyr
#' @import magrittr
#' @import data.table
#'
#' @export

oncoKB_tables <- function(table_names="drug_associations"){
  out_list  <- list()
  tb_nm_opts<- c("drug_associations","therapies","variants","levels","oncogenes","all")
  if(!all(sapply(table_names,function(x) x %in% tb_nm_opts))) stop("Table names can include ",oxford_collapse(tb_nm_opts),".")
  if("all" %in% tolower(table_names)) table_names <- tb_nm_opts
  if("therapies" %in% tolower(table_names)){
    out_list[["therapies"]] <- system.file("extdata","25JAN31-fda_approved_oncology_therapies_uni.csv",package="clugPac") %>%
      fread(col.names = c("year_approval","fda_drug","fda_biomarker","class","mech_or_target","genes","is_targeted","is_precision_onco","ngs_compatible_biomarker")) %>%
    as_tibble %>%
    mutate_at(c("is_targeted","is_precision_onco","ngs_compatible_biomarker"), function(x){
      x=="Y"
    }) %>%
    select(fda_drug,fda_biomarker,class,mech_or_target,genes,is_targeted,is_precision_onco,year_approval,everything())
  }
  if("drug_associations" %in% tolower(table_names)){
    out_list[["drug_associations"]] <- system.file("extdata","25JAN30-oncokb_biomarker_drug_associations.tsv",package="clugPac") %>%
      fread(col.names=c("level","gene_name","alterations","cancer_types","drugs")) %>%
      as_tibble
  }
  if("variants" %in% tolower(table_names)){
    out_list[["variants"]] <- system.file("extdata","25JAN30-oncokb_fda_recognized_alterations.tsv",package="clugPac") %>%
      fread(col.names=c("level","gene_name","alterations","cancer_types")) %>%
      as_tibble
  }
  if("levels" %in% tolower(table_names)){
    out_list[["levels"]] <- system.file("extdata","25JAN30-oncokb_level_ids.txt",package="clugPac") %>%
      fread(col.names=c("level","definition")) %>%
      as_tibble
  }
  if("oncogenes" %in% tolower(table_names)){
    tb_out  <- system.file("extdata","25JAN30-oncokb_cancer_gene_list.tsv",package="clugPac") %>%
      fread(header=TRUE) %>%
      as_tibble %>%
      dplyr::rename_all(function(x) gsub("[ -]","_",gsub("#","num",gsub("[()]","",tolower(x))))) %>%
      dplyr::rename(gene_name = hugo_symbol)
    bool_cols <- tb_out[1,] %>% unlist %>% sapply(function(x) x %in% c("Yes","No")) %>% which
    out_list[["oncogenes"]] <- mutate_at(tb_out,bool_cols, function(x) x == "Yes")
  }
  if(length(out_list) == 1){
    out_list <- out_list[[1]]
  }
  return(out_list)
}

# tb_oncokb_all <- oncoKB_tables("all")
# tb_oncokb_all$drug_associations %>%
#   filter(!grepl("[\\+ ]",drugs)) %>%
#   select(drugs) %>%
#   unlist() %>%
#   unique() %>%
#   sort
#
# tb_oncokb_all$therapies %>%
#   rowwise %>%
#   mutate(genes = str_split(genes,pattern=";")) %>%
#   unnest(genes) %>%
#   filter(genes != "") %>%
#   filter(!genes %in% gr_gns$gene_name) %>%
#   select(genes) %>% unlist %>% unique
#
# gr_gns[grepl("RET",gr_gns$gene_name)]$gene_name %>% sort
# gr_gns[gr_gns$gene_id == "ENSG00000066405"]
