#' @title 
#' Gene biotype color palette.
#' 
#' @description 
#' Returns the ggsci::pal_igv() color palette for all gene biotypes.
#' 
#' @import ggsci
#'
#' @export

biotype_palette   <- function(){
  bts   <- c('IG_C_gene','IG_C_pseudogene','IG_D_gene',
             'IG_J_gene','IG_J_pseudogene','IG_pseudogene',
             'IG_V_gene','IG_V_pseudogene','lncRNA','miRNA',
             'misc_RNA','polymorphic_pseudogene','processed_pseudogene',
             'protein_coding','pseudogene','ribozyme','rRNA',
             'rRNA_pseudogene','scaRNA','scRNA','snoRNA','snRNA',
             'sRNA','TEC','TR_C_gene','TR_D_gene','TR_J_gene',
             'TR_J_pseudogene','TR_V_gene','TR_V_pseudogene',
             'transcribed_processed_pseudogene','transcribed_unitary_pseudogene',
             'transcribed_unprocessed_pseudogene','translated_processed_pseudogene',
             'translated_unprocessed_pseudogene','unitary_pseudogene',
             'unprocessed_pseudogene','vault_RNA')
  pal   <- setNames(ggsci::pal_igv()(length(bts)),bts)
  return(pal)
}