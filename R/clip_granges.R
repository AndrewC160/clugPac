#' @title
#' Clip GRange
#'
#' @description
#' Trim features within a genomic range to fit within a specified window,
#' clipping feature starts/ends to fit within that region. If <replace_cols> is
#' TRUE, start and end columns are changed. Otherwise, start_adj and end_adj
#' columns are appended as metadata and the input GRange object is otherwise
#' unchanged. Can also be provided a tibble which will be converted to a GRange
#' assuming GenomicRanges::makeGRangesFromDataFrame() can interpret the columns
#' present.
#'
#' @param grange_in GenomicRange object(s) to be truncated. Can be a tibble, in which case it will be converted.
#' @param grange_window GenomicRange representing the window to be clipped to. Must be a single range.
#' @param return_as_tibble Should results be returned as a tibble? Defaults to FALSE.
#' @param replace_cols Should start and end values be changed to clipped values? Defaults to FALSE.
#' @param include_clipped_col Should a column be appended to metadata with the number of bases clipped? Defaults to FALSE.
#'
#' @import GenomicRanges
#' @import IRanges
#' @import tidyr
#' @import tibble
#' @import magrittr
#' @import dplyr
#'
#' @export

clip_granges  <- function(grange_in,grange_window,return_as_tibble=TRUE,replace_cols = FALSE,include_clipped_col=FALSE){
  mutate  <- dplyr::mutate
  if(length(grange_window) > 1) stop("subset_grange expects a single GRange for <grange_window>.")
  if(inherits(grange_in,"data.frame")) grange_in  <- makeGRangesFromDataFrame(grange_in,keep.extra.columns = TRUE)

  x_rng   <- c(start=BiocGenerics::start(grange_window),end=BiocGenerics::end(grange_window))

  tb_out  <- subsetByOverlaps(grange_in,grange_window,ignore.strand=TRUE) %>%
    as_tibble %>%
    mutate(
      start_adj = ifelse(start < x_rng['start'],x_rng['start'],start),
      end_adj = ifelse(end > x_rng['end'],x_rng['end'],end))
  if(include_clipped_col)  tb_out <- mutate(tb_out,start_clipped=start_adj - start,end_clipped=end - end_adj)
  if(replace_cols) tb_out <- mutate(tb_out,start = start_adj,end = end_adj) %>% select(-start_adj,-end_adj)
  if(!return_as_tibble)  tb_out <- makeGRangesFromDataFrame(tb_out,keep.extra.columns = TRUE)
  return(tb_out)
}

