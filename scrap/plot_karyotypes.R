#' @title
#' Plot karyotypes
#'
#' @description
#' Plot karyotyping facet for a given genomic range.
#'
#' @param grange_in GRanges object to illustrate.
#'
#' @import ggplot2
#' @import GenomicRanges
#' @import dplyr
#' @import tibble
#' @import magrittr
#'
#' @export

plot_karyotypes <- function(grange_in=NULL){
  if(is.null(grange_in)){
    if(length(grange_in) > 1){
      message("More than one region provided for karyogram, only first GRange will be retrieved.")
      grange_in <- gr[1]
    }
  }
  tb_seq_adj <- enframe(clugPac::get_seqsizes(),name = "seqnames",value="size") %>%
    dplyr::mutate(
      size = as.double(size),
      adj_start = cumsum(lag(size,default = 0))) %>%
    dplyr::mutate(shade = ifelse(row_number() %% 2,"clear","shade"))
  seq_adj <- vectify(tb_seq_adj,adj_start,seqnames)

  tb_bands  <- get_karyotypes(grange_in) %>%
    keepStandardChromosomes(pruning.mode = "coarse") %>%
    as_tibble %>%
    dplyr::filter(seqnames != "chrM") %>%
    dplyr::mutate(
      adj_start = start + seq_adj[as.character(seqnames)],
      adj_end = adj_start + width)
  x_rng   <- c(min(tb_bands$adj_start),max(tb_bands$adj_end))

  ggplot(tb_bands,aes(xmin=adj_start,xmax=adj_end,ymin=0,ymax=1,fill=fill)) +
    scale_x_continuous(limits=x_rng,expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_identity() +
    geom_rect(color=NA) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text.y = element_text(angle=0,hjust=0,vjust=0.5),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size=0.5,color="black",fill=NA),
          plot.margin = unit(c(0,0,0,0),units = "lines"))
}
