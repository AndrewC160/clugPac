#' @title
#' Get seqsizes and adjust to end-to-end x-coordinates.
#'
#' @description
#' Retrieve chromosome sizes for hg38 genome, then calculate adjusted X values
#' that allow chromosomes to be placed end-to-end (for instance, the "adjusted"
#' start of chr1 is 0, that of chr2 is the length of 1, that of chr3 is the sum
#' of the lengths of 1 and 2, and so on).
#'
#' @param seqsize_file Filename of alternate seqsize TSV file to be used. Defaults to NULL, in which case hg38_seqsizes.tsv is loaded from extdata.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#'
#' @export

get_seqsizes_adj  <- function(seqsize_file=NULL){
  get_seqsizes(seqsize_file) %>%
    enframe(name="seqnames",value="size") %>%
    mutate(size=as.double(size),
           adj_start = dplyr::lag(cumsum(size),default = 0)) %>%
    vectify(adj_start,seqnames) %>%
    return()
}
