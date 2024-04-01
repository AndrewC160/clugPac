#' @title BigWig track
#'
#' @description
#' Given a BigWig file and a region of interest, load BigWig data into a
#' table format for plotting.
#'
#' @param bw_file BigWig file of interest.
#' @param gr_window Region of interest in GRanges format.
#' @param bw_colnames Column names to retrieve in BigWig file; defaults to "score".
#' @param bin_num Number of bins to break data into, if desired. Defaults to 1,000, but can also be NULL for no binning.
#'
#' @import rtracklayer
#' @import data.table
#' @import dplyr
#' @import tidyr
#'
#' @export

read_bw_track <- function(bw_file,gr_window,bw_colnames="score",new_colname=NULL,bin_num=1000){
  first <- dplyr::first
  rename<- dplyr::rename

  if(is.null(new_colname)) new_colname  <- bw_colnames

  bws <- BigWigSelection(ranges = gr_window,colnames=bw_colnames)
  bwd <- import.bw(bw_file,selection = bws,as="GRanges") %>%
    as_tibble %>%
    rename(score = !!as.name(bw_colnames)) %>%
    select(-width,-strand) %>%
    mutate(center = (start + end) / 2)

  if(!is.null(bin_num)){
    bwd <- bwd %>%
      select(-start,-end) %>%
      mutate(bin = cut(center,breaks=bin_num,include.lowest=TRUE,dig.lab=12)) %>%
      group_by(bin) %>%
      summarize(seqnames=first(seqnames),
                score = mean(score),
                .groups="drop") %>%
      mutate(start = str_match(as.character(bin),"^.([[:digit:]\\.-]+)")[,2] %>% as.double %>% round(digits=0),
             end = str_match(as.character(bin),",([[:digit:]\\.-]+)")[,2] %>% as.double %>% round(digits=0),
             center = (start + end) / 2) %>%
      select(seqnames,start,end,score,center)
  }
  bwd <- rename(bwd,!!as.name(new_colname) := score)
  return(bwd)
}
