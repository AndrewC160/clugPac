#' @title Localize OGM Reads
#' 
#' @description
#' To visualize OGM contigs relative to their local alignments, pseudo-alignment
#' values can be used to illustrate where the rest of the contig falls were it 
#' to be placed in a given genomic region. The localize_ogm_reads() function 
#' filters an OGM table to only include contigs with at least one XMAP alignment
#' to a genomic region of interest, and also appends columns to the output table
#' which include local pseudo-alignments. For instance, if a 1MB contig consists
#' of a locally-aligned 500MB XMAP and another 500MB XMAP aligned to a different
#' chromosome, the local XMAP column values will mirror its actual coordinates,
#' but the other XMAP will be assigned values that allow it to be drawn in the 
#' same location. In the event that multiple XMAPs align locally, a contig is 
#' returned for each along with a warning (unless 'quiet'=TRUE). Handling these 
#' may need to be done on a case-by-case basis. 
#' 
#' Added columns include:
#' 
#' local: T/F; is the tag part of a locally-aligned XMAP?
#' 
#' contig_xmap: Unique ID comprising contig_id_xmap_id.
#' 
#' contig_strand`: Direction of the entire contig *once aligned to the local 
#' XMAP*. In other words, if the local XMAP is in the minus direction the entire
#' contig is aligned in the minus direction.
#' 
#' local_contig_start: Local (pseudo) genomic position at which the full contig 
#' starts. Will be larger than local_contig_end if the contig is aligned in the 
#' 3'-5' direction.
#' 
#' local_contig_end: Local (pseudo) genomic position at which the full contig 
#' ends.
#' 
#' local_xmap_start: Local (pseudo) genomic position at which a tag should be 
#' drawn. Should equal the local xmap_start/end values for locally aligned tags.
#' 
#' @param gr_window GRanges describing the window to be searched.
#' @param table_in OGM table (as produced by localize_ogm_reads()) to be filtered and formatted.
#' @param quiet Boolean; should warnings about contigs with multiple locally aligned XMAPs be shown? Defaults to TRUE.
#' 
#' @import tidyr
#' @import magrittr
#' @import GenomicRanges
#' @import dplyr
#' @import tibble
#' 
#' @export

localize_ogm_reads <- function(gr_window,table_in,quiet=FALSE){
  rename <- dplyr::rename
  mutate <- dplyr::mutate
  arrange<- dplyr::arrange
  filter <- dplyr::filter
  
  # Get contigs with at least one XMAP local to the gr_window.
  tb_contigs <- filter(table_in,!is.na(start)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE,end.field="start") %>%
    subsetByOverlaps(gr_window) %>%
    as_tibble %>%
    select(contig_id,xmap_id) %>%
    distinct
  
  # Check for and warn about contigs with multiple locally-aligned XMAPs.
  tb_multi_align <- tb_contigs %>% 
    group_by(contig_id) %>%
    mutate(n = n()) %>%
    filter(n > 1)
  
  if(nrow(tb_multi_align) > 0 & !quiet){
    message("Multiple local XMAPs in contigs ",paste(unique(as.character(tb_multi_align$contig_id)),collapse = ", "),".\nMultiple contigs returned.")
  }
  
  # Filter table to keep all contigs with at least one local XMAP.
  tb_x <- table_in %>%
    filter(contig_id %in% tb_contigs$contig_id) %>%
    mutate(local = xmap_id %in% tb_contigs$xmap_id) %>%
    select(local,everything()) %>%
    arrange(contig_id,site_id)
  
  # Get local pseudo-alignment coordinates (start/ends of contigs) for each 
  # contig/XMAP pair (cases where multiple XMAPs align locally are calculated
  # independently).
  tb_loc <- tb_x %>% 
    filter(local) %>%
    select(contig_id,xmap_id,strand,xmap_start,xmap_end,contig_start,contig_end,contig_len) %>%
    distinct %>%
    mutate(local_contig_start = 
             ifelse(strand == "+",
                    xmap_start - contig_start,
                    xmap_end + contig_end),
           local_contig_end = 
             ifelse(strand == "+",
                    local_contig_start + contig_len,
                    local_contig_start - contig_len)) %>%
    rename(contig_strand = strand,
           local_xmap = ) %>%
    mutate(contig_xmap = paste(contig_id,xmap_id,sep="_")) %>%
    select(contig_xmap,contig_id,contig_xmap,contig_strand,local_contig_start,local_contig_end)
  
  tb_x <- tb_x %>%
    left_join(tb_loc,by="contig_id") %>%
    mutate(local_xmap_start = ifelse(contig_strand == "+",
                                     local_contig_start + contig_start,
                                     local_contig_start - contig_start),
           local_xmap_end = ifelse(contig_strand == "+",
                                   local_contig_start + contig_end,
                                   local_contig_start - contig_end),
           local_tag_start = ifelse(contig_strand == "+",
                                    local_contig_start + tag_contig_start,
                                    local_contig_start - tag_contig_start))
  return(tb_x)
}