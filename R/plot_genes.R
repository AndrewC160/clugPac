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
#' @param label_size Size of text in labels' font.
#' @param label_seed Seed number to use for ggrepel labels (used to maintain consistent label / text positions). Defaults to NULL.
#' @param text_genes Gene names to be labeled with text. Defaults to NULL, in which case no genes are labeled.
#' @param text_size Size of text font.
#' @param text_biotypes Vector of gene biotypes to be labeled with text. Defaults to "protein_coding" and "lncRNA." If "all" is listed, all biotypes will be labeled.
#' @param text_seed Seed number to use for ggrepel text (used to maintain consistent label / text positions). Defaults to NULL.
#'
#' @import ggplot2
#' @import magrittr
#' @import GenomicRanges
#' @import dplyr
#' @import ggrepel
#'
#' @export

plot_genes<- function(genomic_region,gr_genes=NULL,label_genes=NULL,label_size=3,label_seed=NULL,text_genes=NULL,text_size=2,text_biotypes=c("protein_coding","lncRNA"),text_seed=NULL){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  x_rng   <- c(start=GenomicRanges::start(genomic_region),end=GenomicRanges::end(genomic_region))
  if(is.null(gr_genes)){
    message("Gathering genes from GFF file, provide a GR object to 'gr_genes' in order to speed this up in the future.")
    gr_genes <- gtf_to_genes() %>% keepStandardChromosomes(pruning.mode="coarse")
  }
  gr <- gr_genes[queryHits(findOverlaps(gr_genes,genomic_region))]

  base_plot   <-
    ggplot(data=NULL) +
    scale_x_continuous(name = grange_desc(genomic_region),
                       limits=x_rng,expand=c(0,0),labels = comma) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(linewidth=0.5,color="black",fill=NA),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linewidth=0.25,color="gray",linetype = "dotted"),
          axis.text.x = element_text(face="italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank())

  if(length(gr) > 0){
    if(is.null(label_genes))label_genes <- ""
    if(is.null(text_genes)) text_genes  <- ""
    if("all" %in% tolower(text_biotypes)) text_biotypes  <- unique(gr$gene_biotype)
    tb_p  <- gr %>%
      as_tibble %>%
      mutate(ymin = pile_coords(start_vals = start,end_vals = end,seqname_vals = seqnames),
             ymax = ymin + 1,
             start_adj = ifelse(start < x_rng[1],x_rng[1],start),
             end_adj = ifelse(end > x_rng[2],x_rng[2],end),
             start = ifelse(strand == "-",end_adj,start_adj),
             end = ifelse(strand == "-",start_adj,end_adj)) %>%
      select(-start_adj,-end_adj) %>%
      mutate(start2 = start + 0.2 * (end - start),
             end2 = end - 0.2 * (end - start),
             lab_type = case_when(gene_name == "" ~ "none",
                                  gene_name %in% label_genes ~ "label",
                                  gene_name %in% text_genes ~ "text",
                                  gene_biotype %in% text_biotypes ~ "text",
                                  TRUE ~ "none"))
    y_rng <- c(0,max(c(6,max(tb_p$ymax))))

    base_plot  <- base_plot +
      scale_y_continuous(limits=y_rng) +
      geom_rect(data=tb_p,
                mapping=aes(xmin = start,xmax=end,
                            ymin = ymin, ymax=ymax,
                            fill = gene_biotype),
                size=0.25,color="black",alpha=0.7) +
      geom_segment(data=filter(tb_p,strand != "*"),
                   mapping=aes(x=start2,xend = end2,
                               y=(ymin + ymax)/2,yend=(ymin + ymax) / 2),
                   size=1,color="black",
                   arrow = arrow(length = unit(0.5,"lines"),type = "closed"))

    #Gene text.
    tb_labs <- filter(tb_p,lab_type == "text")
    if(nrow(tb_labs) > 0){
      base_plot <- base_plot +
        geom_text_repel(data=tb_labs,
                  seed = text_seed,
                  size=text_size,
                  mapping=aes(x=(start + end) / 2,
                              y=(ymin + ymax) / 2,
                              label = gene_name),
                  min.segment.length = 0.01,
                  nudge_y = sample(c(1,-1),replace=TRUE))#,size = nrow(tb_labs)))
    }
    #Gene labels.
    tb_labs <- filter(tb_p,lab_type == "label")
    if(nrow(tb_labs) > 0){
      base_plot <- base_plot +
        geom_label_repel(
          data=tb_labs,
          seed = label_seed,
          mapping=aes(x=(start + end) / 2,
                      y=(ymin + ymax) / 2,
                      label = gene_name),
          min.segment.length = 0.01,
          size=label_size,
          nudge_y = 2)
    }
  }
  return(base_plot)
}
