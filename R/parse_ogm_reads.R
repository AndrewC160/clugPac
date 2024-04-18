#' @title Parse OGM reads
#'
#' @description
#' OGM data can be read from a combination of XMAP and CMAP files, which contain
#' contigs and their alignments to the genome called by Bionano and the
#' locations of fluorescent tags on each contig, respectively. The
#' `parse_ogm_reads()` function allows this to be done with pairs of XMAP/CMAP
#' files, and includes a name column which either contains the input XMAP
#' filename or, if multiple samples are read at the same time and at least one
#' of the two lists is named, the names assigned to the input XMAP/CMAP file
#' lists. This function reads pairs of files and compiles them into a single
#' table representing contigs (the actual DNA structure assembled from
#' overlapping reads), XMAPs (alignments assigned to each region of a given
#' contig), and tags (fluorescent tag locations). Since tags are only defined by
#' their position in a contig, this function calculates their putative *genomic*
#' location by determining which XMAP they fall in and measuring their position
#' based on its coordinates. Note that while XMAPs may overlap (there can be
#' ambiguity for alignments in some contigs), tag genomic locations are
#' calculated after assigning them to the highest-confidence XMAP. Also note
#' that tag locations may differ slightly as XMAP alignments are confidence
#' based.
#'
#' The output table uses one row per tag with each tag's calculated genomic
#' location as the `seqnames` and `start` positions. Other fields are as follows:
#'
#' strand: XMAP's purported alignment direction.
#'
#' contig_id/xmap_id/site_id: ID values for the entire contig, each reported
#' genomic alignment *within* that contig (can be overlapping), and each
#' fluorescent tag location, respectively.
#'
#' coverage: Total reads that comprise the contig which include each tag. I.e. a
#' coverage of 3 indicates that 3 of the molecules assembled into that contig
#' included this tag, so higher coverage suggests greater confidence.
#'
#' contig_len: Length of the total contig in bp.
#'
#' contig_start/end: Start/end positions of each XMAP alignment *within the
#' contig* in bp. I.e., an XMAP with contig_start of 2,000 starts 2kb from the
#' start of the contig. Note that contig_start is greater than contig_end in
#' minus-strand XMAPs.
#'
#' xmap_start/end: *Genomic* start/end positions for each XMAP according to
#' Bionano.
#'
#' xmap_conf: Confidence score of XMAP according to Bionano.
#'
#' tag_contig_start: Location of tag *within the contig* in bp, i.e. a tag with
#' tag_contig_start of 2,000 is located 2kb from the start of the contig.
#'
#' name: Name of the sample in question; defaults to the XMAP file name if only
#' one pair of XMAP/CMAP files are provided, but if a named list of files are
#' provided these names will be used.
#'
#' @param fls_xmap Character vector or list of XMAP files to read. Can be named, and its names will take precedence over names of fls_cmap.
#' @param fls_cmap Character vector or list of CMAP files to read. Can be named.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import magrittr
#' @import GenomicRanges
#'
#' @export

parse_ogm_reads <- function(fls_xmap,fls_cmap){
  rename <- dplyr::rename
  mutate <- dplyr::mutate
  arrange<- dplyr::arrange
  filter <- dplyr::filter

  if(length(fls_xmap) > 1){
    if(length(fls_xmap) != length(fls_cmap)) stop(paste0("Different number of XMAP/CMAP files provided (",length(fls_xmap)," vs. ",length(fls_cmap),", respectively)."))
    if(is.null(names(fls_xmap))){
      if(is.null(names(fls_cmap))){
        names(fls_xmap) <- fls_xmap
      }
      names(fls_xmap) <- names(fls_cmap)
    }
  tb_xmap <- lapply(1:length(fls_xmap), function(i){
    nm <- names(fls_xmap)[i]
    parse_ogm_reads(fls_xmap =fls_xmap[i],
                    fls_cmap =fls_cmap[i]) %>%
      mutate(name = nm) %>%
      select(name,everything())
    }) %>%
      do.call(rbind,.)
  }else{
    # Read XMAP alignments.
    tb_xmap <- read_xmap_file(fls_xmap,as_granges = FALSE) %>%
      mutate(contig_id = as.character(contig_id),
             xmap_id = as.character(xmap_id))

    # Read CMAP alignments.
    tb_cmap <- read_cmap_file(fls_cmap) %>%
      mutate(contig_id = as.character(contig_id)) %>%
      mutate(tag_id = paste0(contig_id,"_",site_id))

    # Assign tags to XMAPs: Convert XMAP and CMAP files into GenomicRanges objects
    # using contig_ids as seqnames, use annotate_gr() to annotate tags with their
    # assigned XMAPs. Note that only the HIGHEST CONFIDENCE XMAP is assigned.
    gr_xmap <- tb_xmap %>%
      arrange(desc(conf)) %>%
      mutate(seqnames = contig_id,
             start = ifelse(contig_start > contig_end,contig_end,contig_start),
             end = ifelse(contig_start > contig_end,contig_start,contig_end)) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    gr_cmap <- tb_cmap %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                               seqnames.field = "contig_id",
                               start.field = "position",
                               end.field = "position")

    tb_cmap <- annotate_gr(gr_cmap,gr_xmap,cols_query = "xmap_id",cols_subject="xmap_id",first_only = TRUE) %>%
      as_tibble %>%
      rename(contig_id = seqnames,
             tag_contig_start = start) %>%
      select(contig_id,xmap_id,site_id,tag_contig_start,coverage)

    tb_xmap <- inner_join(tb_xmap,tb_cmap,by=c("contig_id","xmap_id")) %>%
      mutate(contig_id = factor(contig_id)) %>%
      mutate(xmap_len = contig_end - contig_start) %>%
      mutate(tag_xmap_frac = (tag_contig_start - contig_start)/xmap_len) %>%
      mutate(tag_xmap_frac = ifelse(strand == "-",1-tag_xmap_frac,tag_xmap_frac)) %>%
      mutate(tag_xmap_pos = tag_xmap_frac * xmap_len) %>%
      mutate(tag_genome = ifelse(strand == "+",
                                 start + tag_xmap_pos,
                                 start + tag_xmap_pos - xmap_len)) %>%
      rename(xmap_start = start,
             xmap_end = end,
             start = tag_genome,
             xmap_conf = conf) %>%
      select(seqnames,start,strand,contig_id,xmap_id,site_id,coverage,contig_len,
             contig_start,contig_end,xmap_start,xmap_end,
             xmap_conf,tag_contig_start) %>%
      arrange(seqnames,start) %>%
      mutate(name=fls_xmap)
  }
  return(tb_xmap)
}
