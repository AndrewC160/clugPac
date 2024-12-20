#' @title Plot genomic RNA
#'
#' @description
#' Given a genomic range of gene locations (output by gtf_to_genes(), for
#' instance), and columns containing sample names, gene names, and TPM values,
#' plot RNA expression as scatter plots (if fewer than <max_points> genes are
#' present in the region) or tiles. Also outputs aligned gene tracks, by default.
#'
#' @param gr_window_in GenomicRanges region to plot.
#' @param gr_genes_in GRanges of gene locations.
#' @param sample_names_in Vector of sample names for each TPM value.
#' @param sample_genes_in Vector of gene names for each TPM value.
#' @param table_tpm_in Vector of TPM values.
#' @param plot_metric Value to plot. Defaults to TPM, can also be "ratio" in which case the ratio of expression among samples is plotted.
#' @param z_scale Scale for values. Defaults to log2, but can also be "log10."
#' @param max_points Maximum number of genes to display as scatter plots; defaults to 20.
#' @param point_size Point size for scatter plots.
#' @param include_gene_track Type of gene track to include. Defaults to "genes", but can also be "pileup" or "none."
#'
#' @import dplyr
#' @import tidyr
#' @import gridExtra
#' @import grid
#' @import ggplot2
#' @import GenomicRanges
#' @import tibble
#'
#' @export

plot_rna_genome <- function(gr_window_in,gr_genes_in,sample_names_in,sample_genes_in,sample_tpm_in,
                            plot_metric="TPM",z_scale="log2",max_points=20,point_size=2,include_gene_track="genes",
                            widths_in=c(1,12,1),heights_in=c(1,5,5,1),top_trace=TRUE,bottom_trace=FALSE){
  tb_rna <- tibble(gene=sample_genes_in,
                   sample_name=sample_names_in,
                   tpm=sample_tpm_in) %>%
    mutate(sample_name = factor(sample_name,levels=select(.,sample_name) %>% unlist %>% unique))

  tb_p <- subsetByOverlaps(gr_genes_in,gr_window_in) %>%
    as_tibble %>%
    inner_join(tb_rna,by=c("gene_name"="gene"),relationship="many-to-many") %>%
    mutate(gene_pos = (end + start)/2) %>%
    arrange(gene_pos) %>%
    mutate(gene_name = factor(gene_name,levels=select(.,gene_name) %>% unlist %>% unique)) %>%
    group_by(gene_name) %>%
    mutate(x_idx = as.integer(gene_name),
           gene_pos_frac = (gene_pos-start(gr_window_in)) / width(gr_window_in)) %>%
    ungroup %>%
    group_by(seqnames,start,end,gene_name,sample_name,x_idx,gene_pos_frac) %>%
    summarize(tpm = mean(tpm),.groups="drop") %>%
    mutate(y_idx = as.integer(sample_name)) %>%
    ungroup %>%
    group_by(gene_name) %>%
    mutate(ratio = tpm / median(tpm,na.rm=TRUE)) %>%
    rowwise %>%
    mutate(x_jitter = x_idx + runif(n=1,min = -0.25,max=0.25)) %>%
    ungroup

  if(tolower(plot_metric) == "ratio"){
    tb_p  <- mutate(tb_p,score=ratio)
    z_lab <- "TPM ratio"
  }else if(tolower(plot_metric) == "tpm"){
    tb_p  <- mutate(tb_p,score=tpm)
    z_lab <- "TPM"
  }else{
    stop("<plot_metric> must be either 'TPM' or 'Ratio'.")
  }
  if(z_scale=="log2"){
    tb_p  <- mutate(tb_p,score=log2(score))
    z_lab <- paste0("log2(",z_lab,")")
  }else if(z_scale == "log10"){
    tb_p  <- mutate(tb_p,score=log10(score))
    z_lab <- paste0("log10(",z_lab,")")
  }

  x_rng <- c(0.5,max(tb_p$x_idx)+0.5)

  tb_segs <- tb_p %>%
    group_by(gene_name) %>%
    summarize(x_idx=dplyr::first(x_idx),
              mean_tpm = mean(tpm),
              gene_pos_frac = dplyr::first(gene_pos_frac),
              .groups='drop') %>%
    mutate(gene_pos = x_rng[1] + gene_pos_frac * diff(x_rng))

  x_labs <- tb_p %>%
    group_by(gene_name) %>%
    summarize(x_idx = dplyr::first(x_idx)) %>%
    mutate(shade = x_idx %% 2 == 0)

  if(length(unique(tb_p$gene_name)) > max_points){
    y_rng <- c(-0.5,max(tb_p$y_idx)+0.5)
    trace_segs_top <- NULL
    if(top_trace){
      y_rng <- y_rng + c(0,0.5)
      trace_segs_top <-
        geom_segment(data=tb_segs,mapping=aes(x=x_idx,xend=gene_pos),
                     y=max(tb_p$y_idx)+0.5,yend=y_rng[2],inherit.aes=FALSE,
                     color="gray85",show.legend = FALSE)
    }
    if(bottom_trace){
      trace_bottom_color <- "gray85"
    }else{
      trace_bottom_color <- "white"
    }
    trace_segs_bot <-
      geom_segment(data=tb_segs,mapping=aes(x=x_idx,xend=gene_pos),
                   y=0.5,yend=y_rng[1],inherit.aes=FALSE,
                   color=trace_bottom_color,show.legend = FALSE)
    y_labs<- tb_p %>%
      group_by(sample_name) %>%
      summarize(y_idx = dplyr::first(y_idx)) %>%
      vectify(y_idx,sample_name)

    p <- ggplot(tb_p,aes(xmin=x_idx-0.5,xmax=x_idx+0.5,
                    ymin=y_idx-0.5,ymax=y_idx+0.5,
                    fill=score,label=gene_name)) +
      scale_x_continuous(limits=x_rng,expand=c(0,0)) +
      scale_y_continuous(name="RNA",breaks=as.double(y_labs),labels=names(y_labs),
                         expand=c(0,0),limits=y_rng) +
      scale_fill_viridis(name=z_lab) +
      geom_rect() +
      trace_segs_top +
      trace_segs_bot +
      geom_text(data=x_labs,mapping=aes(x=x_idx,label=gene_name),y=0.5,inherit.aes=FALSE,
                angle=90,hjust=1,vjust=0.5,size=2) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_blank(),
            plot.margin=unit(c(0,0,0,0),"lines"),
            panel.background = element_blank(),
            plot.background = element_rect(fill="white",color=NA),
            panel.grid = element_blank())
  }else{
    y_rng <- c(0,max(tb_p$score)*1.1,8)
    trace_segs_top <- NULL
    trace_segs_bot <- NULL
    if(top_trace){
      y_rng <- y_rng + c(0,0.1)
      trace_segs_top <-
        geom_segment(data=tb_segs,mapping=aes(x=x_idx,xend=gene_pos),
                   y=max(y_rng),yend=max(y_rng)*1.1,inherit.aes=FALSE,
                   color="gray85",show.legend = FALSE)
    }
    if(bottom_trace){
      y_rng <- y_rng + c(-0.1,0)
      trace_segs_bot <-
        geom_segment(data=tb_segs,mapping=aes(x=x_idx,xend=gene_pos),
                     y=min(y_rng),yend=min(y_rng)*1.1,inherit.aes=FALSE,
                     color="gray85",show.legend = FALSE)
    }
    #y_rng <- c(0,max(1.1 *max(tb_p$score),8))
    y_brks<- ceiling(seq(0,ceiling(max(tb_p$score)),length.out=4))

    p <- ggplot(tb_p,aes(x=x_jitter,y=score,color=sample_name,label=gene_name)) +
      scale_x_continuous(breaks=x_labs$x_idx,labels=x_labs$gene_name,expand=c(0,0),
                         limits=x_rng,oob = scales::oob_keep) +
      scale_y_continuous(name=z_lab,expand=c(0,0),limits=1.1*y_rng,breaks=y_brks) +
      scale_alpha_manual(values=c(`TRUE`=0.2,`FALSE`=0)) +
      annotate(geom="rect",xmin=x_rng[1],xmax=x_rng[2],ymin=y_rng[1],ymax=y_rng[2],fill="white",color="black",linewidth=0.5) +
      geom_rect(data=x_labs,mapping=aes(xmin=x_idx-0.5,xmax=x_idx+0.5,alpha=shade),
                ymin=y_rng[1],ymax=y_rng[2],fill="blue",color=NA,linewidth=0.5,
                inherit.aes=FALSE,show.legend = FALSE) +
      trace_segs_top +
      trace_segs_bot +
      geom_point(size=point_size) +
      theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=6),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_blank(),
            plot.background = element_rect(fill="white",color=NA),
            plot.margin=unit(c(0,0,0,0),"lines"),
            panel.background = element_blank(),
            legend.key = element_rect(fill=NA,color=NA),
            legend.box = element_blank(),
            legend.title = element_blank(),
            panel.grid = element_blank())
  }
  if(include_gene_track!="none"){
    p_genes <- NULL
    if(include_gene_track=="pileup"){
      p_genes <- plot_genes_pileup(genomic_region = gr_window_in,gr_genes=gr_genes_in)
    }else if(include_gene_track=="genes"){
      p_genes <- plot_genes(genomic_region = gr_window_in,gr_genes = gr_genes_in,
                              text_bins = 40,x_trace_alpha = 0.3,
                              text_genes = as.character(sample_genes_in),
                              x_trace_genes = unique(as.character(tb_p$gene_name)))
    }else{
      stop("<include_gene_track> should be 'none', 'pileup', or 'genes'.")
    }

    p_genes <- plot_split(p_genes +
                            theme(plot.margin = unit(c(0,0,-0.2,0),"lines"),
                                  plot.background = element_rect(fill="white",color=NA),
                                  legend.position = "none",
                                  axis.text.x = element_text(),
                                  axis.title.x = element_text()))
    p <- plot_split(p)

    p <- grid.arrange(
      p_genes$bottom,
      p_genes$left,p_genes$main,p_genes$right,
      p$left,p$main,p$right,p$bottom,
      layout_matrix=
        rbind(c(NA,1,NA),
              c(2,3,4),
              c(5,6,7),
              c(NA,8,NA)),
      widths=widths_in,
      heights=heights_in)
  }
  return(p)
}
