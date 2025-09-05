#' @title Plot transcripts
#'
#' @description
#' Given a contiguous genomic region of interest and a GTF file, plot gene
#' transcripts and transcript features from that GTF file.
#'
#' @param genomic_region Congiguous genomic region to illustrate.
#' @param gtf_file GTF file of feature annotations. If NULL, default GTF from gtf_extract_features() is used.
#' @param focus_genes Gene names/Ensembl IDs to highlight. Defaults to NULL, in which case all genes are given an alpha of 1.
#' @param gene_text_size Size of text on gene annotations. Defaults to 3.
#' @param exon_text_size Size of text on exon annotations. Defaults to 2.
#' @param max_support_level Maximum transcript support level (1 being highest confidence). Defaults to 7, which is any transcript support level.
#' @param max_transcripts Maximum number of transcripts to display for a given gene even if <max_support_level> filter would includ more. Defaults to 7, in which case the 7 highest confidence transcripts are shown.
#' @param label_exons Should exons be labeled with their exon number? Defaults to FALSE.
#' @param return_tables Should plotting be skipped and formatted tables be returned instead? Defaults to FALSE.
#'
#' @import dplyr
#' @import tidyr
#' @import GenomicRanges
#' @import ggplot2
#' @import scales
#' @import ggrepel
#' @import magrittr
#'
#' @export

plot_transcripts  <- function(genomic_region,gtf_file=NULL,focus_genes=NULL,gene_text_size=3,exon_text_size=2,max_support_level = 7,max_transcripts=7,label_exons=FALSE,return_tables=FALSE){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  arrange <- dplyr::arrange

  x_rng   <- c(start=GenomicRanges::start(genomic_region),end=GenomicRanges::end(genomic_region))

  if(length(genomic_region) > 1) stop("plot_transcripts() is intended for contiguous regions only.")

  gr_feats<- gtf_extract_features(genomic_region,gtf_file = gtf_file,retain_attributes=c("gene_name","gene_id","transcript_id","exon_id","exon_number","transcript_support_level","gene_version"))
  seqlevelsStyle(gr_feats) <- "UCSC"
  gr_feats<- clugPac::clip_granges(grange_in = gr_feats,grange_window = genomic_region,include_clipped_col = TRUE)

  focus_gns <- focus_genes %||% c(unique(gr_feats$gene_name),unique(gr_feats$gene_id))

  scale_x <- scale_x_continuous(name = grange_desc(genomic_region),
                                oob = scales::oob_keep,
                                limits=x_rng,expand=c(0,0),labels = comma)

  tb_sum  <- gr_feats %>%
    as_tibble %>%
    filter(!is.na(transcript_id)) %>%
    rename(support=transcript_support_level) %>%
    mutate(support=ifelse(support %in% as.character(1:5),support,"6"),
           support=as.integer(support)) %>%
    group_by(seqnames,gene_name,gene_id,transcript_id) %>%
    summarize(gene_min = min(start),
              gene_max=max(end),
              t_scripts = length(unique(transcript_id)),
              support = min(na.omit(support)),
              width= sum(width),
              .groups='drop') %>%
    filter(support <= max_support_level) %>%
    arrange(gene_name,support,desc(width),gene_min,gene_max) %>%
    group_by(gene_name,gene_id) %>%
    filter(row_number() <= max_transcripts) %>%
    ungroup() %>%
    mutate(y = pile_coords(gene_min,gene_max),
           focus = gene_name %in% focus_gns | gene_id %in% focus_gns) %>%
    select(gene_name,gene_id,transcript_id,gene_min,gene_max,support,focus,y)

  tb_p <- gr_feats %>%
    as_tibble %>%
    select(-width,-score,-frame) %>%
    inner_join(tb_sum,by=c("gene_name","gene_id","transcript_id")) %>%
    mutate(feature = factor(feature,levels = c("transcript",
                                               "exon",
                                               "CDS",
                                               "three_prime_utr",
                                               "five_prime_utr",
                                               "stop_codon",
                                               "start_codon"))) %>%
    filter(!is.na(feature))

  if(return_tables){
    p <- tb_p
  }else{
    feat_cols <- c(transcript="gray70",
                   exon = "black",
                   CDS = "blue",
                   three_prime_utr="red",
                   five_prime_utr ="green",
                   stop_codon = "darkred",
                   start_codon="darkgreen")

    geom_ex_labs <- NULL
    if(label_exons){
      geom_ex_labs <- geom_text_repel(
        data=filter(tb_p,focus & feature == "exon"),
        mapping=aes(label=exon_number),
        min.segment.length = 0.001,size = exon_text_size)
    }

    p <- ggplot(tb_p,aes(xmin=start,xmax=end,
                    x=start,xend=end,
                    ymin=y+0.2,ymax=y+0.8,
                    y=y+0.5,yend=y+0.5,
                    label=gene_name,
                    fill=feature,
                    color=feature)) +
      scale_x +
      scale_color_manual(name="Features",values=feat_cols) +
      scale_fill_manual(name="Features",values=feat_cols) +
      # facet_wrap(.~support,ncol=1,scales="free_y",strip.position = "right") +
      # scale_alpha_manual(name="Support",values=c("1"=1,
      #                                            "2"=0.8,
      #                                            "3"=0.5,
      #                                            "4"=0.3,
      #                                            "5"=0.2,
      #                                            "6"=0.5)) +
      geom_segment(data=filter(tb_p,feature == "transcript"),show.legend=FALSE) +
      geom_rect(data=filter(tb_p,!feature %in% c("transcript","stop_codon","start_codon")),
                color=NA) +
      geom_text(data=filter(tb_p,feature == "transcript" & strand == "-"),hjust=1,show.legend=FALSE,size=gene_text_size) +
      geom_text(data=filter(tb_p,feature == "transcript" & strand == "+"),hjust=0,show.legend=FALSE,size=gene_text_size,
                mapping=aes(x=end)) +
      geom_point(data=filter(tb_p,feature %in% c("start_codon","stop_codon")),pch=4) +
      geom_ex_labs +
      theme(plot.background = element_rect(fill="white",color=NA),
            panel.border = element_rect(fill=NA,color="black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle=45,hjust=1,vjust=1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank())
  }
  return(p)
}
