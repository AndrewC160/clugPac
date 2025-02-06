#' @title Binned CNA average
#'
#' @description
#' Given a genomic range object of interest and another genomic range object
#' representing genomic bins with a numeric column for copy number values,
#' calculate the average copy number of the input ranges by averaging all bins
#' overlapped with respect to total area covered. For instance, consider a
#' feature 10kb wide which overlaps 10 1kb bins. If the coverage value of 9 bins
#' is 1 and the remaining bin is 10, the average coverage of the feature becomes
#' 1.9:
#'
#' 90% of the feature's area is 1 and 10% is 10, so 0.9 x 1 + 0.1 x 10 = 1.9.
#'
#' The input range will be returned with an extra column containing the coverage
#' values. The column of <gr_cna_in> with coverage values is specified using
#' <cna_col>, and if an output column name is not provided via <query_col>, this
#' will be the column name used for the output. Note that if an existing column
#' is found of the same name, the function stops with an error.
#'
#' @param gr_in GRanges object to be quantified.
#' @param gr_cna_in GRanges object of genome-wide bins containing coverage data.
#' @param cna_col Column name in <gr_cna_in> containing coverage values.
#' @param query_col Column name to be given to coverage values in the output GRanges object; defaults to the value of <cna_col>.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import GenomicRanges
#' @import GenomicAlignments
#'
#' @export

binned_cna_average <- function(gr_in,gr_cna_in,cna_col,query_col){
  if(missing(query_col)) query_col <- cna_col
  if(query_col %in% colnames(as_tibble(gr_in))) stop("'",query_col,"' column already present in <gr_in>.")

  gr_cna<- gr_cna_in[,NULL]
  olaps <- findOverlaps(gr_in,gr_cna)
  cvg   <- width(
    pintersect(gr_in[queryHits(olaps)],
               gr_cna[subjectHits(olaps)]))
  wid   <- width(gr_in[queryHits(olaps)])
  frac  <- cvg / wid
  cna   <- mcols(gr_cna_in)[subjectHits(olaps),cna_col]
  tot   <- frac * cna

  tb_cna<- olaps %>%
    as_tibble %>%
    cbind(cna=tot) %>%
    group_by(queryHits) %>%
    summarize(cna = sum(cna),.groups='drop')

  gr_out<- gr_in %>%
    as_tibble %>%
    mutate(queryHits = row_number()) %>%
    left_join(tb_cna,by="queryHits") %>%
    select(-queryHits) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  if(any(is.na(gr_out$cna))) warning("NAs detected in copy number column; were some regions not covered by CNA data?")

  return(gr_out)
}
