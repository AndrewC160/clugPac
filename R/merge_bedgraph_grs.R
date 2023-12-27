#' @title Merge BedGraph GRs
#' 
#' @description
#' Given two GRanges with bedGraph information, combine them and average their
#' score columns. Function uses pintersect() so the output range will typically
#' include additional features. For instance:
#' 
#' grA:   ====4====|====2======
#' grB:   =2=|========4========
#' grOut: =3=|==4==|=====3=====
#' 
#' @param grA_in GRanges A.
#' @param grB_in GRanges B.
#' @param score_column Score column name. Defaults to "score".
#' 
#' @import GenomicRanges 
#' 
#' @export

merge_bedgraph_grs <- function(grA_in,grB_in,score_column="score"){
  grA <- grA_in[,score_column]
  grB <- grB_in[,score_column]
  olaps <- findOverlaps(grA,grB)
  grAB <- pintersect(grA[queryHits(olaps)],
                     grB[subjectHits(olaps)])
  scrsA <- mcols(grA)[[score_column]][queryHits(olaps)]
  scrsB <- mcols(grB)[[score_column]][subjectHits(olaps)]
  mcols(grAB)[[score_column]] <-(scrsA + scrsB)/2
  gr_out <- grAB[width(grAB) > 1]
  return(gr_out)
}
