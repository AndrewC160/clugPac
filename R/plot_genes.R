#' @title
#' Plot gene track
#'
#' @description
#' Given a genomic range, plots a simple gene track using gene annotations from
#' gtf_to_genes(). Returns a GGplot2 object.
#'
#' @param genomic_region GRange denoting the region to be illustrated.
#' @param gr_genes GRanges object containing gene annotations. If not provided, gtf_to_genes() is used, but this can be time consuming: better to use gtf_to_genes() earlier and provide the output GR to this function directly.
#' @param label_genes Gene names to be labeled. Defaults to NULL, in which case no genes are labeled. Supersedes text_genes.
#' @param label_size Size of text in labels' font. Defaults to 3.
#' @param label_seed Seed number to use for ggrepel labels (used to maintain consistent label positions). Defaults to NULL.
#' @param focus_genes Gene names to highlight (i.e. give an alpha of 1). If given, genes not in this list have their alpha reduced to <background_alpha>. Defaults to NULL, in which case all genes are alpha 1.
#' @param background_alpha Alpha value to use for background (non-<focus-genes>). Defaults to 0.25.
#' @param text_genes Gene names to be labeled with text. Defaults to NULL, in which case no genes are labeled. Genes must have non-empty gene name columns.
#' @param text_biotypes Vector of gene biotypes to be labeled with text. Defaults to "protein_coding" and "lncRNA." If "all" is listed, all biotypes will be labeled.
#' @param text_size Size of text font. Defaults to 2.
#' @param text_bins X-axis is split into <text_bins> "bins" and labels are stacked within each. Fewer bins result in higher stacks of labels, while lower bin numbers may result in more overlaps with labels/geometry. Defaults to 10.
#' @param text_nudge_y Distance above each text bin's highest gene annotation above which to start adding text labels in units of gene layers (one layer = one gene annotation). Defaults to 0.25.
#' @param text_line_height Height of text labels, useful to adjust in cases where larger text sizes are used causing them to overlap. Units of layers, defaults to 0.2.
#' @param text_nudge_x Distance to nudge labels in the negative x-direction when staggering stacked labels. Units of fractions of the X-axis, defaults to 0.005 (text is nudged 0.5% of the x-axis per stacked label).
#'
#' @import ggplot2
#' @import magrittr
#' @import GenomicRanges
#' @import dplyr
#' @import ggrepel
#' @import scales
#' @import ggsci
#'
#' @export

# gr_genes <- gtf_to_genes() %>% keepStandardChromosomes(pruning.mode="coarse")
# genomic_region <- gr_genes[gr_genes$gene_name == "DLX5"] %>% resize(width=5E6,fix="center")
# text_biotypes="all"
# label_genes <- "DLX5"

plot_genes<- function(genomic_region,gr_genes=NULL,label_genes=NULL,label_size=3,label_seed=NULL,
                      focus_genes=NULL,background_alpha=0.25,
                      text_genes=NULL,text_biotypes=c("protein_coding","lncRNA"),
                      text_size=2,text_bins=10,text_nudge_y=0.25,text_line_height=0.2,text_nudge_x = 0.005){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  x_rng   <- c(start=GenomicRanges::start(genomic_region),end=GenomicRanges::end(genomic_region))
  if(is.null(gr_genes)){
    message("Gathering genes from GFF file, provide a GR object to 'gr_genes' in order to speed this up in the future.")
    gr_genes <- gtf_to_genes() %>% keepStandardChromosomes(pruning.mode="coarse")
  }
  gr <- gr_genes[queryHits(findOverlaps(gr_genes,genomic_region))]

  base_plot   <- ggplot(data=NULL) +
    scale_x_continuous(name = grange_desc(genomic_region),
                       limits=x_rng,expand=c(0,0),labels = comma) +
    scale_fill_igv(name="Biotype") +
    guides(alpha="none",
           color="none") +
    theme(plot.background = element_rect(fill="white",color=NA),
          plot.margin = unit(c(0,0,0,0),"lines"),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth=0.5,color="black",fill=NA),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linewidth=0.5,color="gray",linetype = "dotted"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())

  if(length(gr) > 0){
    if(is.null(label_genes))label_genes <- ""
    if(is.null(text_genes)) text_genes  <- ""
    if(is.null(focus_genes)) focus_genes<- unique(gr$gene_name)
    if("all" %in% tolower(text_biotypes)) text_biotypes  <- unique(gr$gene_biotype)
    bt_levs <- c("protein_coding","miRNA","lncRNA")
    bt_levs <- c(bt_levs,setdiff(levels(gr_genes$gene_biotype),bt_levs))
    tb_p  <- gr %>%
      as_tibble %>%
      mutate(ymin = pile_coords(start_vals = start,end_vals = end,seqname_vals = seqnames,min_gap=0.001 * diff(x_rng)),
             ymax = ymin + 1,
             start_adj = ifelse(start < x_rng[1],x_rng[1],start),
             end_adj = ifelse(end > x_rng[2],x_rng[2],end),
             start = ifelse(strand == "-",end_adj,start_adj),
             end = ifelse(strand == "-",start_adj,end_adj),
             gene_biotype = factor(gene_biotype,levels=bt_levs)) %>%
      select(-start_adj,-end_adj) %>%
      mutate(start2 = start + 0.2 * (end - start),
             end2 = end - 0.2 * (end - start),
             mid = (start + end) / 2,
             focus_type = ifelse(gene_name %in% focus_genes,"foreground","background"),
             lab_type = case_when(gene_name == "" ~ "none",
                                  gene_name %in% label_genes ~ "label",
                                  gene_name %in% text_genes ~ "text",
                                  gene_biotype %in% text_biotypes ~ "text",
                                  TRUE ~ "none"))

    tb_labs <- filter(tb_p,lab_type != "none")
    if(nrow(tb_labs) > 0){
      tb_labs<- tb_labs %>%
        arrange(desc(start),desc(end)) %>%
        mutate(text_bin = cut(start,breaks=text_bins,dig.lab = 10)) %>%
        mutate(text_bin_start = str_match(as.character(text_bin),"^[\\[\\(]([[:digit:]\\.]+)")[,2] %>% as.double %>% ceiling,
               text_bin_end = str_match(as.character(text_bin),"([[:digit:]\\.]+)[\\]\\)]$")[,2] %>% as.double %>% ceiling) %>%
        group_by(text_bin) %>%
        mutate(text_y1 = max(ymax) + text_nudge_y,
               text_y2 = text_y1 + text_line_height* row_number() * max(ymax),
               text_bin_start = text_bin_start - (row_number() * text_nudge_x * diff(x_rng)),
               text_bin_start = ifelse(text_bin_start < x_rng[1],x_rng[1],text_bin_start))
    }else{
      tb_labs <- mutate(tb_labs,text_y1 = 0,text_y2 = 0)
    }

    # Scale y up for very large windows to represent "zoomed out" scale.
    min_y_scale <- 6 + ceiling(width(genomic_region)/1E6)/2
    y_rng <- c(0,max(c(min_y_scale,tb_labs$text_y2,tb_p$ymax)))
    scale_y <- scale_y_continuous(name = "Genes",limits=c(-1,y_rng[2]+1),expand=c(0,0))

    #Gene text.
    tb_l <- filter(tb_labs,lab_type == "text")
    if(nrow(tb_l) > 0){
      base_plot <- base_plot +
        geom_segment(data=tb_l,
                     mapping=aes(x=text_bin_start,
                                 xend=text_bin_start,
                                 y=text_y2,yend=text_y1,
                                 alpha=focus_type),
                     color="gray75") +
        geom_segment(data=tb_l,
                     mapping=aes(x=text_bin_start,
                                 xend=start,
                                 y=text_y1,yend=ymax,
                                 alpha=focus_type),
                     color="gray75") +
        geom_segment(data=tb_l,
                     mapping=aes(x=start,
                                 xend=start,
                                 y=ymax,yend=(ymin+ymax)/2,
                                 alpha=focus_type),
                     color="gray75") +
        geom_text(data=tb_l,
                  mapping=aes(x=text_bin_start,
                              y=text_y2,
                              label=gene_name,
                              alpha=focus_type),
                  hjust=0,vjust=0,
                  size=text_size)
    }

    base_plot  <- base_plot +
      scale_alpha_manual(values=c(foreground=1,background=background_alpha)) +
      scale_color_manual(values=c(foreground="black",background=NA)) +
      scale_y +
      geom_rect(data=tb_p,
                mapping=aes(xmin = start,xmax=end,
                            ymin = ymin+0.1, ymax=ymax-0.1,
                            fill = gene_biotype,
                            alpha=focus_type,
                            color=focus_type),
                linewidth=0.25) +
      geom_segment(data=filter(tb_p,strand != "*" & focus_type != "background" & lab_type != "none"),
                   mapping=aes(x=start2,xend = end2,
                               y=(ymin + ymax)/2,yend=(ymin + ymax) / 2,
                               alpha=focus_type,
                               color=focus_type),
                   linewidth=1,
                   arrow = arrow(length = unit(0.5,"lines"),type = "closed"))

    #Gene labels.
    tb_l <- filter(tb_labs,lab_type == "label")
    if(nrow(tb_l) > 0){
      base_plot <- base_plot +
        geom_label_repel(
          data=tb_l,
          seed = label_seed,
          mapping=aes(x=(start + end) / 2,
                      y=ymin+0.1,
                      label = gene_name,
                      alpha=focus_type),
          min.segment.length = 0.01,
          size=label_size,ylim = c(0,-0.5))
    }
  }
  return(base_plot)
}
