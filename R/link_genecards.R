#' @title Link to GeneCards
#' 
#' @description Link a given gene name to a GeneCards search: Basically pastes
#' the gene name/ID into a GeneCards URL and hopes for the best. Does not check
#' if links work, so should be reasonable when generating columns of a table.
#' 
#' @param gene_name Name or ID of gene to link to. Should work with ENS IDs, etc.
#' @param link_to_search Should link direct to a search page instead of the actual gene page? Defaults to FALSE.
#' 
#' @export

link_genecards  <- function(gene_name = "six2",link_to_search = FALSE){
  if(link_to_search){
    #https://www.genecards.org/Search/Keyword?queryString=SIX2
    lnk_txt <- paste0("https://www.genecards.org/Search/Keyword?queryString=",tolower(gene_name))
  }else{
    #https://www.genecards.org/cgi-bin/carddisp.pl?gene=SIX2&keywords=SIX2
    lnk_txt <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene_name)
  }
  return(lnk_txt)
}
