#' @title Chromosome color palette
#' 
#' @description
#' Returns a named vector of colors from ggsci::ucscgb useful for consistent
#' chromosome colors.
#' 
#' @import ggsci
#' 
#' @export

chrom_colors <- function(){
  setNames(ggsci::pal_ucscgb()(24),paste0("chr",c(1:22,"X","Y")))
}
