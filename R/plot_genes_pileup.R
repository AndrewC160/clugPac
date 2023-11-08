#' @title
#' Plot gene track pileup.
#'
#' @description
#' Given a genomic range, breaks the range into <bin_num> bins and plots the
#' number of genes per bin. Specific genes can be labeled using <label_genes>.
#'
#' @param genomic_region GRange denoting the region to be illustrated.
#' @param gr_genes_in GRanges object containing gene annotations. If not provided, gtf_to_genes() is used, but this can be time consuming: better to use gtf_to_genes() earlier and provide the output GR to this function directly.
#' @param bin_num Number of bins to split the range into. Defaults to 300.
#' @param label_genes Gene names to be labeled. Defaults to NULL, in which case no genes are labeled. Supersedes text_genes.
#' @param label_size Size of text in labels' font. Defaults to 3.
#' @param label_seed Seed number to use for ggrepel labels (used to maintain consistent label positions). Defaults to NULL.
#'
#' @import ggplot2
#' @import magrittr
#' @import GenomicRanges
#' @import dplyr
#' @import ggrepel
#' @import scales
#'
#' @export

plot_genes_pileup <- function(genomic_region,gr_genes_in=NULL,bin_num=300,label_genes=NULL,label_size=3,label_seed=NULL,kary_frac=0.1){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  x_rng   <- c(start=GenomicRanges::start(genomic_region),end=GenomicRanges::end(genomic_region))
  if(is.null(gr_genes_in)){
    message("Gathering genes from GFF file via gtf_to_genes(), provide a GR object to 'gr_genes_in' in order to speed this up in the future.")
    gr_genes_in <- gtf_to_genes() %>% keepStandardChromosomes(pruning.mode="coarse")
  }
  gr_gns <- subsetByOverlaps(gr_genes_in,genomic_region)

  gr_bins <- tile(genomic_region,n = bin_num) %>% unlist
  gr_bins$bin_id  <- paste0("bin_",c(1:length(gr_bins)))
  gr_bins$count <- countOverlaps(gr_bins,gr_genes_in)

  tb_p <- as_tibble(gr_bins) %>%
    mutate(pMb = count / (width/1e6))
  bin_size <- prettyBP(mean(width(gr_bins)),digits=0)

  y_max <- max(tb_p$count)
  y_rng <- NULL
  geom_labs <- NULL
  lab_frac  <- 0.2 + kary_frac
  if(!is.null(label_genes)){
    tb_labs <- gr_gns[gr_gns$gene_name %in% label_genes] %>% as_tibble
    if(nrow(tb_labs) > 0){
      y_rng <- c(-y_max*lab_frac,y_max*1.1)
      geom_labs <- geom_label_repel(data=tb_labs,mapping=aes(y=-kary_frac * y_max,x=(start+end)/2,label=gene_name),
                                    size=label_size,nudge_y = -kary_frac * y_max,inherit.aes=FALSE,
                                    min.segment.length = 0)
    }
  }
  if(is.null(y_rng)){
    y_rng <- c(-y_max * kary_frac,y_max*1.1)
  }

  tb_kary <- get_karyotypes() %>%
    subsetByOverlaps(genomic_region) %>%
    as_tibble %>%
    mutate(ymax=0,ymin=-y_max*kary_frac)

  ggplot(tb_p,aes(xmin=start,xmax=end,ymin=0,ymax=count)) +
    scale_x_continuous(name=grange_desc(genomic_region,append_ending = paste0("; ",bin_size," bin size")),
                       expand=c(0,0),limits=x_rng,labels=prettyBP) +
    scale_y_continuous(name="Genes per bin",expand=c(0,0),limits=y_rng) +
    geom_rect(fill="gray40",color="gray20") +
    geom_rect(data=tb_kary,mapping=aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax,fill=fill),
              inherit.aes=FALSE) +
    annotate(geom="rect",xmin=x_rng[1],xmax=x_rng[2],ymax=0,ymin=-kary_frac*y_max,fill=NA,color="black",linewidth=0.5) +
    scale_fill_identity() +
    geom_labs +
    theme(plot.background = element_rect(fill="white",color=NA),
          plot.margin = unit(c(0,0,0,0),"lines"),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth=0.5,color="black",fill=NA),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linewidth=0.5,color="gray",linetype = "dotted"))
}
