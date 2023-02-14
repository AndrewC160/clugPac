#' @title
#' Downsample BDG table.
#'
#' @description
#' Given a table containing BedGraph format data, downsample genomic
#' coordinates such that approximately <bin_num> bins are included. Useful for
#' plotting (etc.) when BedGraphs of several megabases are to be shown on an
#' axis that will ultimately only span a few thousand pixels, both in terms of
#' storage and processing time. Binning is done separately for each chromosome
#' by splitting the the table by chromosome and running function on each
#' subset. Defaults to mean scores for each bin, but any meaningful function
#' can be provided via <bin_function>, "bin_function = max" for instance.
#'
#' @param table_in Table of BedGraph values. Must have standard start/end/seqnames/score columns.
#' @param bin_num Number of bins to split coords into. Defaults to 2,000. If NULL or shorter than the actual range, the original table will be returned.
#' @param name_col Column containing names. Defaults to "name".
#' @param bin_function Function to combine bin values. Defaults to mean, but can be any function that takes in multiple values and returns a single number (sum, max,etc.).
#'
#' @import dplyr
#' @import stringr
#' @import tidyr
#'
#' @export

downsample_bdg_table  <- function(table_in,bin_num = 2000,name_col="name",bin_function=mean){
  if(is.null(bin_num)){
    return(table_in)
  }
  if(length(unique(table_in$seqnames)) > 1){
    #If more than one chromosome is included, do each separately.
    tb_out <- table_in %>%
      as_tibble %>%
      dplyr::group_split(seqnames) %>%
      lapply(downsample_bdg_table) %>%
      do.call(rbind,.)
    return(tb_out)
  }
  if(!is.null(name_col) & !is.null(table_in[,name_col])){
    tb_in <- rename(table_in,name = !!as.name(name_col))
  }else{
    tb_in <- mutate(table_in,name = "")
  }
  x_rng   <- c(min(table_in$start),max(table_in$end))

  if(diff(x_rng) <= bin_num){
    table_out   <- table_in
  }else{
    brks  <- seq(x_rng[1],x_rng[2],length.out = bin_num)
    names(brks) <- c(1:bin_num)

    tb_out<- tb_in %>%
      mutate(mid = (start + end) / 2,
             bin = cut(mid,breaks = brks,include.lowest = TRUE,dig.lab = 12)) %>%
      group_by(seqnames,name,bin) %>%
      summarize(score = bin_function(score),.groups="drop") %>%
      mutate(start=str_match(as.character(bin),"^[\\[\\(]([:digit:]+)")[,2] %>% as.integer,
             end = str_match(as.character(bin),"([[:digit:]\\.]+)[\\)\\]]$")[,2] %>% as.integer) %>%
      select(seqnames,start,end,score,name)

    if(!is.null(name_col)){
      tb_out  <- rename(tb_out,!!as.name(name_col) := name)
    }else{
      tb_out  <- select(tb_out,-name)
    }
  }
  return(tb_out)
}
