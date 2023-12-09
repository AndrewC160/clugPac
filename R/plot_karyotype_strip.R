#' @title
#' Plot karyotype track.
#'
#' @description
#' Plot genome-wide or chromosome-specific karytotypes. Plot labels/title are
#' included, but hidden by default using theme(). If <window_region> is set to
#' "chrom" or "genome", the strip will include the entire chromosome or genome
#' (respectively) and the <genomic_region> will be illustrated with a red box
#' and two red arrows.
#'
#' @param genomic_region GRange denoting the region to be illustrated.
#' @param window_region Region to be illustrated overall; defualts to "chrom" (plots entire chromosome), and can be "genome" (plots entire genome). Other values cause only the selected region to be drawn.
#' @param trace_line_color Line color for box to draw illustrating the <genomic_region>.
#' @param trace_fill Fill of trace; defaults to "red".
#' @param trace_point_size Size of arrows on trace; defaults to 7.
#' @param trace_line_width Linewidth of trace; defaults to 0.5.
#' @param trace_alpha Alpha of trace; defaults to zero (empty red rectangle).
#' @param trace_linetype Linetype of trace; defaults to "solid".
#' @param line_color_in Color of lines around genomic bands; defaults to NA (no lines).
#' @param line_widths_in Width of lines around genomic bands; defaults to 0.5.
#' @param label_size_in Size of labels to draw; defaults to 5.
#' @param label_frac_in Fraction of x-axis below which bands are not labeled; defaults to 0.02, in which case bands which occupy less than 2% of the x-axis are drawn but not labeled.
#'
#' @import ggplot2
#' @import magrittr
#' @import GenomicRanges
#' @import dplyr
#' @import scales
#' @import ggnewscale
#'
#' @export

# genomic_region <- GRanges("chr15",IRanges(42124999,43124999))

plot_karyotype_strip  <- function(genomic_region,window_region="chrom",trace_line_color="red",trace_point_size=7,trace_line_width=0.5,trace_alpha=0,trace_fill="red",trace_linetype="solid",line_color_in=NA,line_widths_in=0.5,label_size_in=5,label_frac_in=0.02){
  start <- BiocGenerics::start
  end   <- BiocGenerics::end
  gr_kary <- get_karyotypes()
  gr_seqs <- get_seqsizes(as_granges=TRUE)
  gr_chrom<- gr_seqs[as.character(seqnames(gr_seqs)) == as.character(seqnames(genomic_region))]
  gr_adj  <- get_seqsizes_adj()[as.character(seqnames(gr_chrom))]

  y_rng   <- c(0,1)
  y_pad   <-0.01
  tb_trace<- as_tibble(genomic_region) %>%
    mutate(start_adj = start + gr_adj,
           end_adj = end + gr_adj,
           ymin=y_rng[1],ymax=y_rng[2],
           label="")
  scale_x <- NULL
  if(window_region == "genome"){
    x_rng   <- c(1,max(gr_kary$end_adj))
    tb_kary <- as_tibble(gr_kary) %>%
      mutate(label="")
    x_ttl   <- "Genome"
    scale_x <- scale_x_continuous(name=x_ttl,
                                  breaks=mcols(gr_seqs)$adj_start,
                                  labels=as.character(seqnames(gr_seqs)),
                                  expand=c(0,0),
                                  oob = scales::oob_squish)
  }else if(window_region == "chrom"){
    x_rng   <- c(start(gr_chrom),end(gr_chrom)) + gr_adj
    tb_kary <- clip_granges(gr_kary,gr_chrom) %>%
      mutate(start_adj = start + gr_adj,
             end_adj = end + gr_adj)
    x_ttl   <- grange_desc(gr_chrom)
  }else{
    x_rng   <- c(start(genomic_region),end(genomic_region)) + gr_adj
    tb_kary <- clip_granges(gr_kary,genomic_region) %>%
      mutate(start_adj = start + gr_adj,
             end_adj = end + gr_adj)
    x_ttl   <- grange_desc(genomic_region)
    tb_trace<- NULL
  }
  if(is.null(scale_x)){
    scale_x <- scale_x_continuous(name=x_ttl,limits=x_rng,expand=c(0,0),labels=prettyBP,
                                  oob = scales::oob_squish)
  }
  tb_kary   <- tb_kary %>%
    mutate(ymin=y_rng[1]+y_pad,ymax=y_rng[2]-y_pad,
           x_frac = width / sum(width)) %>%
    mutate(label = ifelse(x_frac > label_frac_in,name,""),
           alpha=1)

  p <- ggplot(tb_kary,mapping=aes(xmin=start_adj,xmax=end_adj,
                             ymin=ymin,ymax=ymax,
                             x = (start_adj + end_adj)/2,
                             y = diff(y_rng)/2,
                             fill=fill,color=color,
                             label=label)) +
    scale_x +
    scale_y_continuous(limits=y_rng,expand=c(0,0)) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_linetype_identity() +
    scale_linewidth_identity() +
    scale_alpha_identity() +
    geom_rect(color=line_color_in,linewidth=line_widths_in) +
    geom_text(angle=45,hjust=0.5,vjust=0.5) +
    annotate(geom="rect",xmin=x_rng[1],xmax=x_rng[2],ymin=y_rng[1]+y_pad,ymax=y_rng[2]-y_pad,
             fill=NA,color="black") +
    gg_theme("genome") +
    ggtitle(x_ttl) +
    theme(plot.background = element_rect(fill="white",color=NA),
          plot.margin = unit(c(0,0,0,0),"lines"),
          plot.title = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  if(!is.null(tb_trace)){
    trace_width <- tb_trace$end[1] - tb_trace$start[1]
    reg_width   <- diff(x_rng)
    if(trace_width <= 0.01 * reg_width){
      p <- p +
        geom_rect(data=tb_trace,color=trace_line_color,fill=trace_fill,alpha=trace_alpha,linetype=trace_linetype,linewidth=trace_line_width)
    }
    p <- p +
      geom_point(data=tb_trace,mapping = aes(x=(start_adj+end_adj)/2,y=ymax-y_pad/2),
                 pch = 25,color="black",alpha=1,size=trace_point_size,fill=trace_fill) +
      geom_point(data=tb_trace,mapping = aes(x=(start_adj+end_adj)/2,y=ymin+y_pad/2),
                 pch = 24,color="black",alpha=1,size=trace_point_size,fill=trace_fill)
  }
  return(p)
}
