#' @title
#' Get karyotypes.
#'
#' @description
#' Get karyotype data from stored cytobands.txt file. Retreives either the full
#' genome or all bands within a specified genomic range. Includes "suggested"
#' colors in the "fill" column.
#'
#' @param grange_in GRanges object denoting area to retrieve. Assumes a contiguous region; for non-contiguous regions use multiple calls.
#' @param standard_chroms Should only standard chromosomes be included? Defaults to TRUE.
#'
#' @import GenomicRanges
#' @import dplyr
#' @import tibble
#' @import magrittr
#' @import data.table
#'
#' @export

get_karyotypes  <- function(grange_win=NULL,standard_chroms=TRUE){
  fill_table <- c(gneg="white",
                   gpos25="gray75",
                   gpos50="gray50",
                   gpos75="gray25",
                   gpos100="black",
                   acen = "red",
                   stalk = "lightblue",
                   gvar = "pink")
  color_table<- c(gneg="black",
                  gpos25="black",
                  gpos50="black",
                  gpos75="white",
                  gpos100="white",
                  acen = "black",
                  stalk = "black",
                  gvar = "black")
  sq_adj  <- get_seqsizes_adj()
  gr <- fread(system.file("extdata","cytobands.txt",package = "clugPac"),sep="\t") %>%
    dplyr::rename(
      seqnames=`#chrom`,
      start = chromStart,
      end = chromEnd) %>%
    as_tibble %>%
    dplyr::mutate(start = ifelse(start ==0,1,start),
                  fill = fill_table[gieStain],
                  color= color_table[gieStain],
                  start_adj = start + sq_adj[as.character(seqnames)],
                  end_adj = end + sq_adj[as.character(seqnames)])
  if(standard_chroms){
    gr <- gr %>%
      mutate(seqnames=factor(seqnames,levels=paste0("chr",c(1:22,"X","Y")))) %>%
      dplyr::filter(!is.na(seqnames))
  }
  gr <- makeGRangesFromDataFrame(gr,seqnames.field = "seqnames",start.field = "start",end.field = "end",keep.extra.columns = TRUE)

  if(!is.null(grange_win)){
    g_min   <- min(start(grange_win) + sq_adj[as.character(seqnames(grange_win))])
    g_max   <- max(end(grange_win) + sq_adj[as.character(seqnames(grange_win))])
    gr  <- as_tibble(gr) %>%
      filter(start_adj <= g_max & end_adj >= g_min) %>%
      mutate(start_off = start_adj - g_min,
             end_off = g_max - end_adj,
             start = ifelse(start_off < 0,start - start_off,start),
             end = ifelse(end_off < 0,end + end_off,end),
             start_adj = ifelse(start_off < 0,start_adj + start_off,start_adj),
             end_adj = ifelse(end_off < 0,end_adj + end_off, end_adj)) %>%
      select(-end_off,-start_off) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  }
  return(gr)
}
