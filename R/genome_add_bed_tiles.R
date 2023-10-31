#' @title
#' Genome: add BED tiles
#'
#' @description
#' Add a layer of tiles (rectangles) with an infinite y-axis range representing
#' start-end coordinates on the genome as stored in a BED file. Named files will
#' be separated into facets, while unnamed files will stretch over all facets.
#'
#' @param plot_base_in Output stemming from <genome_plot()>.
#' @param bed_files List of BED files containing regions to highlight. Named BED files will be printed one-per-facet, unnamed files will be printed across all facets.
#' @param bed_fills List of fill values for ranges. Defaults to "black."
#' @param bed_colors List of color values for ranges. Defaults to "black."
#' @param bed_alphas List of alpha values for ranges. Defaults to 0.3.
#' @param bed_linetypes List of linetypes to apply to ranges. Defaults to "solid."
#'
#' @import GenomicRanges
#' @import IRanges
#' @import tidyr
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#'
#' @export

genome_add_bed_tiles<- function(plot_base_in,bed_files,bed_fills="black",bed_colors="black",bed_alphas=0.3,bed_linetypes="solid",bed_linewidths=0.5){
  tb_fls <-
    tibble(name = names(bed_files),
           fill=bed_fills,
           color=bed_colors,
           alpha=bed_alphas,
           linetype=bed_linetypes,
           linewidth=bed_linewidths,
           bed=bed_files) %>%
    mutate(idx = paste0("bed_",row_number()))

  tb_beds<- apply(tb_fls,1, function(row_in) {
    if("name" %in% colnames(tb_fls)){
      nm <- row_in[['name']]
    }else{
      nm <- "Pileup"
    }
    idx <- row_in[['idx']]
    fread(row_in[['bed']]) %>%
      as_tibble %>%
      dplyr::rename(seqnames=1,start=2,end=3) %>%
      mutate(name = nm,
             idx = idx) %>%
      select(-starts_with("V")) %>%
      return()
  }) %>% do.call(rbind,.) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    clip_granges(grange_window = plot_base_in$gr_win,replace_cols = TRUE,return_as_tibble = TRUE) %>%
    left_join(tb_fls) %>% suppressMessages

  p <- plot_base_in$plot +
    facet_wrap(.~name,ncol=1,scales='free_y',strip.position = "right") +
    geom_rect(data=tb_beds,mapping=aes(fill=fill,color=color,alpha=alpha,linetype=linetype,linewidth=linewidth),
              ymin=-Inf,ymax=Inf)

  plot_base_out <- plot_base_in
  plot_base_out$plot <- p
  return(plot_base_out)
}
