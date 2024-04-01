#' @title Read XMAP file
#' 
#' @description
#' Read a BioNano-formatted XMAP file (file containing OGM 
#' alignment data) into a tibble containing XMAP and contig data. Can also 
#' return a GRanges object.
#' 
#' @param xmap_fl XMAP filename.
#' @param as_grange Boolean; should the table be returned as a GRanges object? Defaults to TRUE.
#' 
#' @import dplyr
#' @import tidyr
#' @import GenomicRanges
#' @import data.table
#' @import stringr
#' @import tibble
#' 
#' @export

read_xmap_file <- function(xmap_fl,as_granges=TRUE){
  rename <- dplyr::rename
  mutate <- dplyr::mutate
  arrange<- dplyr::arrange
  filter <- dplyr::filter
  
  c_nms <- paste("grep '#h'", xmap_fl) %>% 
    system(intern = TRUE) %>% 
    str_split("\t",simplify = TRUE) %>%
    gsub("^#. ","",.)
  
  tb_out <- paste0("grep -v '#' ",xmap_fl) %>%
    fread(cmd = .,col.names=c_nms) %>%
    as_tibble %>%
    mutate(seqnames=RefContigID,
           start=ifelse(RefStartPos > RefEndPos,RefEndPos,RefStartPos),
           end = ifelse(RefStartPos > RefEndPos,RefStartPos,RefEndPos)) %>%
    mutate(seqnames=case_when(seqnames == 23 ~ "X",
                              seqnames == 24 ~ "Y",
                              TRUE ~ as.character(seqnames)) %>%
             paste0("chr",.) %>%
             factor(levels=paste0("chr",c(1:22,"X","Y")))) %>%
    rename(strand=Orientation,
           xmap_id=XmapEntryID,
           contig_id=QryContigID,
           contig_start=QryStartPos,
           contig_end = QryEndPos,
           contig_len = QryLen,
           conf=Confidence) %>%
    select(seqnames,start,end,strand,xmap_id,starts_with("contig"),conf) %>%
    arrange(contig_id,seqnames,start)
  if(as_granges){
    tb_out <- makeGRangesFromDataFrame(tb_out,keep.extra.columns = TRUE)
  }
  return(tb_out)
}