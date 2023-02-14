#' @title
#' Preview GRanges.
#'
#' @description
#' Given a list of GenomicRange objects, this function produces a summary plot
#' using GGplot2, focused on whatever genomic region was provided in
#' "window_coords". If no window_coords provided, attempts to find a region
#' with the most overlaps between all input ranges.
#'
#' @param list_grs List of GenomicRanges, either named or not.
#' @param window_coords GenomicRanges object representing window to focus on. Defaults to NULL.
#'
#' @import GenomicRanges
#' @import IRanges
#' @import ggplot2
#' @import scales
#'
#' @export

preview_granges <- function(list_grs,window_coords=NULL){
  start   <- GenomicRanges::start
  end     <- GenomicRanges::end
  #Convert normal lists to GRangesList.
  list_grs <- GRangesList(list_grs)

  #Each GRange should be named.
  if(is.null(names(list_grs))){
    names(list_grs) <- paste0("Range_",seq_along(list_grs))
  }

  #If one wasn't provided, find a window to view.
  #Randomly select from the first GRanges, prioritize those with most unique
  #overlaps among the subsequent GRanges in the list (if more than one was
  #provided).
  if(is.null(window_coords)){
    if(length(list_grs) > 1){
      #Find regions with at least one overlapping feature among the most GRanges.
      mtx <- sapply(c(2:length(list_grs)), function(i) {
        countOverlaps(list_grs[[1]],list_grs[[i]])
      }) %>%
        as.matrix(nrow=length(list_grs[[1]]))
      mtx <- mtx > 0
      most_laps <- which(rowSums(mtx) == max(rowSums(mtx)))
    }else{
      most_laps <- c(1:length(list_grs[[1]]))
    }
    window_coords   <- list_grs[[1]][sample(most_laps,1)]
    width_val <- width(window_coords)
    if(width_val>1e6){
      width_val <- 1e6
    }else if(width_val < 500){
      width_val <- 1e3
    }else{
      width_val <- 1.1 * width(window_coords)
    }
    window_coords <- resize(window_coords,width = width_val,fix="center")
  }else{
    window_coords <- window_coords[1]
  }
  x_rng <- c(start=start(window_coords),end = end(window_coords))

  #Gather overlapping features
  p_tb <- lapply(c(1:length(list_grs)), function(i) {
    list_grs[[i]][subjectHits(findOverlaps(window_coords,list_grs[[i]]))] %>%
      as_tibble %>%
      select(seqnames,start,end,width,strand) %>%
      mutate(feature = names(list_grs)[i])
  }) %>%
    do.call(rbind,.) %>%
    rowwise %>%
    mutate(feature = factor(feature,levels=names(list_grs)),
           start = max(start,x_rng['start']),
           end = min(end,x_rng['end']))
  x_lab <- grange_desc(grange_in = window_coords)

  ggplot(p_tb,aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=feature)) +
    facet_wrap(.~feature,ncol=1,scales="fixed",strip.position = "right",drop = FALSE) +
    scale_x_continuous(name = x_lab,limits=x_rng,expand=c(0,0),labels=comma) +
    geom_rect(size=0.5,color="black") +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(size=0.5,color="gray",linetype="dotted"),
          panel.spacing.y = unit(0,"cm"),
          strip.background = element_blank(),
          strip.text.y = element_text(angle=0,hjust=0,face="italic"),
          axis.text.x = element_text(face="italic"),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
}
