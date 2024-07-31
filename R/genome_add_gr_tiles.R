#' @title
#' Genome: add BED tiles
#'
#' @description
#' Add a layer of tiles (rectangles) with an infinite y-axis range representing
#' start-end coordinates on the genome as stored in a BED file.
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

genome_add_gr_tiles<- function(plot_base_in,gr_input,gr_name_column=NULL,gr_fills="black",gr_colors="black",gr_alphas=0.3,gr_linetypes="solid",gr_linewidths=0.5){
  if(is.null(gr_name_column)){
    unnamed<- TRUE
    gr_name_column <- "name"
    tb_fls <- tibble(name=gr_name_column,
                     fill=gr_fills[1],
                     color=gr_colors[1],
                     alpha=gr_alphas[1],
                     linetype=gr_linetypes[1],
                     linewidth=gr_linewidths[1],
                     idx=1)
  }else{
    unnamed<- FALSE
    tb_fls <- gr_input %>%
      as_tibble %>%
      mutate(name :=!!as.name(gr_name_column)) %>%
      arrange(name) %>%
      select(name) %>%
      distinct %>%
      mutate(fill=gr_fills,
             color=gr_colors,
             alpha=gr_alphas,
             linetype=gr_linetypes,
             linewidth=gr_linewidths,
             idx = paste0("gr_",row_number()))
  }
  gr_wind <- plot_base_in$gr_win
  tb_regs <- clip_granges(gr_input,grange_window = gr_wind,
                          replace_cols = TRUE,return_as_tibble = TRUE)
  if(unnamed){
    tb_regs <- mutate(tb_regs,name = gr_name_column)
  }
  tb_regs <- tb_regs %>%
    mutate(name := !!as.name(gr_name_column)) %>%
    left_join(tb_fls,by="name")
  if(unnamed) tb_regs <- mutate(tb_regs,name=NULL)
  if(nrow(tb_regs) > 0){
    p <- plot_base_in$plot +
      facet_wrap(.~name,ncol=1,strip.position = "right") +
      geom_rect(data=tb_regs,mapping=aes(fill=fill,color=color,alpha=alpha,linetype=linetype,linewidth=linewidth),
                ymin=-Inf,ymax=Inf)
    plot_base_out <- plot_base_in
    plot_base_out$plot <- p
  }else{
    plot_base_out <- plot_base_in
  }
  return(plot_base_out)
}
