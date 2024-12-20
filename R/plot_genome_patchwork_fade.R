#' @title Plot genome patchwork fade.
#'
#' @description
#' Given a GRanges object denoting a series of genomic regions on the same
#' x-axis (for instance, a contig made up of disparate loci), represent it as a
#' rectangle with fill representing chromosome names. Regions are split into
#' bins and alpha is used to denote strand information.
#'
#' @param gr_regs_in GRanges object of one or more regions to illustrate.
#' @param bin_num Number of bins to split each region into.
#' @param simplify_labels Boolean; should genomic values be written in "10MB" annotation? Defaults to false, i.e. "1,000,000".
#' @param simplify_digits If <simplify_labels>, how many digits to display when rounding. Defaults to 2.
#' @param background_fill Color of the panel background, defaults to "white."
#' @param axis_position Position of the x-axis, defaults to "bottom" but can also be "top".
#'
#' @import tidyr
#' @import magrittr
#' @import tibble
#' @import ggsci
#' @import ggplot2
#'
#' @export

plot_genome_patchwork_fade  <- function(gr_regs_in,bin_num = 40,simplify_labels=TRUE,simplify_digits=2,background_fill="white",axis_position="bottom"){
  select  <- dplyr::select
  filter  <- dplyr::filter
  mutate  <- dplyr::mutate
  rename  <- dplyr::rename

  if(simplify_labels){
    lab_func <- function(x) prettyBP(x,digits=simplify_digits)
  }else{
    lab_func <- comma
  }

  if(is.null(gr_regs_in)){
    gr_regs_in <- get_seqsizes(as_granges=TRUE)
  }

  tb_p <- gr_regs_in %>%
    as_tibble %>%
    select(seqnames,start,end,strand) %>%
    mutate(region = row_number(),
           width=end-start) %>%
    rowwise %>%
    mutate(start2 = list(seq(start,end,length.out=bin_num)),
           x_pos = list(start2 - start),
           idx = list(as.double(1:bin_num))) %>%
    mutate(idx = ifelse(strand == '-',list(rev(idx)),list(idx)),
           x_pos = ifelse(strand == "-",list(rev(x_pos)),list(x_pos))) %>%
    ungroup %>%
    unnest(c(idx,start2,x_pos)) %>%
    group_by(seqnames) %>%
    mutate(seq_min = min(start),
           seq_max = max(end),
           frac = (start2 - start) / (seq_max - seq_min)) %>%
    ungroup %>%
    select(-seq_min,-seq_max) %>%
    arrange(region,idx) %>%
    group_by(region) %>%
    mutate(x_end = lead(x_pos)-1) %>%
    ungroup %>%
    mutate(x_end=ifelse(is.na(x_end),end-start,x_end),
           width= x_end - x_pos,
           x_end = cumsum(width),
           x_start = x_end - width) %>%
    select(seqnames,start,end,start2,strand,region,idx,x_start,x_end,frac)

  tb_labs <- tb_p %>%
    group_by(seqnames,start,end,strand,region) %>%
    summarize(x_start=min(x_start),
              x_end = max(x_end),
              .groups='drop') %>%
    mutate(swap=strand=='-') %>%
    swap_columns("start","end","swap") %>%
    select(-swap) %>%
    mutate(frac=1,
           label = as.character(seqnames),
           label = ifelse(strand == "-",paste0(label,"(-)"),label),
           width = x_end - x_start,
           frac = width / sum(width),
           x_mid = (x_end + x_start) / 2)

  x_ax <- tb_labs %>%
    rowwise %>%
    mutate(x_brks = case_when(frac > 0.7 ~ 5,
                              frac > 0.6 ~ 4,
                              frac > 0.2 ~ 3,
                              TRUE ~ 2)) %>%
    mutate(x_brks = list(seq(x_start,x_end,length.out=x_brks))) %>%
    unnest(x_brks) %>%
    group_by(region) %>%
    mutate(x_frac = (x_brks - min(x_brks)) / (max(x_brks) - min(x_brks)),
           x_vals = (max(end) - min(start))*x_frac + min(start)) %>%
    # In cases where only one pixel is present per region, x-labels break. In
    # these cases manually establish start/end coordinates/fractions.
    mutate(x_frac = case_when(is.na(x_frac) & row_number() == 1 ~ 0,
                              is.na(x_frac) & row_number() == n() ~ 1,
                              TRUE ~ x_frac),
           x_vals = case_when(is.na(x_vals) & row_number() == 1 ~ start,
                              is.na(x_vals) & row_number() == n() ~ end,
                              TRUE ~ x_vals)) %>%
    ##########################################################################
    ungroup %>%
    group_by(region,strand,seqnames) %>%
    summarize(x_brks = list(x_brks),
              x_vals = list(x_vals),
              reg_start = min(start),
              reg_end = max(end),
              .groups="drop") %>%
    rowwise %>%
    mutate(width = reg_end - reg_start,
           x_vals = ifelse(strand == "-",list(rev(x_vals)),list(x_vals))) %>%
    unnest(c(x_brks,x_vals)) %>%
    group_by(region,strand,seqnames) %>%
    mutate(x_frac = (x_vals - min(x_vals)) / (max(x_vals) - min(x_vals)),
           x_vals = reg_start + (width * x_frac),
           bound = case_when(row_number() == 1 ~ "left",
                             row_number() == n() ~ "right",
                             TRUE ~ "mid"),
           label = lab_func(x_vals)) %>%
    ungroup %>%
    mutate(label = ifelse(bound=="right" & lead(bound,default="")=="left",paste(label,lead(label),sep="-\n"),label),
           label = ifelse(bound=="left"& lag(bound,default="")=="right","",label)) %>%
    filter(label != "") %>%
    vectify(x_brks,label)

  ggplot(tb_p,aes(xmin=x_start,xmax=x_end,ymin=0,ymax=1,fill=seqnames,alpha=frac)) +
    scale_x_continuous(expand=c(0,0),breaks=x_ax,labels=names(x_ax),position = axis_position) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=chrom_colors()) +
    geom_rect() +
    geom_rect(data=tb_labs,color="black",fill=NA,linewidth=0.5) +
    geom_text(data=tb_labs,mapping=aes(x=x_mid,y=0.5,label=label),
              alpha=0.9,color='black') +
    theme(plot.background = element_rect(fill="white",color=NA),
          plot.margin = unit(c(0,0,0,0),"lines"),
          panel.background = element_rect(fill=background_fill),
          panel.border = element_rect(fill=NA,color="black",linewidth=0.5),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle=90,vjust=0.5),
          axis.title.x = element_blank(),
          axis.ticks.x = element_line(color="black",linewidth=0.25),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none")
}
