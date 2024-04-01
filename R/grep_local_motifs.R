#' @title Grep Local Motifs
#'
#' @description
#' Given a genomic window of interest, search it for exact matches to a given
#' motif using the 'grep' function (i.e., not a standard motif searching
#' algorithm). Returns a GRanges object of motif locations.
#'
#' @param gr_window GRanges object defining the contiguous region of the genome to search.
#' @param motif Motif in DNA alphabet to search for. Defaults to CTTAAG, which is the DLE1 binding motif for OGM.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import GenomicRanges
#' @import tibble
#' @import Biostrings
#' @import stringr
#'
#' @export

grep_local_motifs <- function(gr_window,motif="CTTAAG"){
  mutate <- dplyr::mutate
  arrange<- dplyr::arrange

  # "Appropriate:" http://localhost:10000/help/library/Biostrings/html/PDict-class.html
  # dict0 <- DNAStringSet(BSgenome.Hsapiens.UCSC.hg38)
  # gen_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38)
  # mtf_seq <- Biostrings::DNAString(motif)
  #
  # p_dct <- PDict(gen_seq,skip.invalid.patterns = TRUE)
  # Biostrings::matchPDict(p_dct,mtf_seq)
  # But I'm fine with exact, p-value-free matches via str_locate_all.
  gen_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr_window) %>% suppressWarnings
  mtf_seq <- Biostrings::DNAString(motif)
  rbind(
    str_locate_all(as.character(gen_seq),
                   as.character(mtf_seq))[[1]] %>%
      as_tibble %>%
      mutate(strand="+"),
    str_locate_all(as.character(gen_seq),
                   as.character(reverseComplement(mtf_seq)))[[1]] %>%
      as_tibble %>%
      mutate(strand="-")) %>%
    mutate(seqnames=as.character(seqnames(gr_window)) %>%
             factor(levels=paste0("chr",c(1:22,"X","Y"))),
           start = start + start(gr_window),
           end = end + start(gr_window)) %>%
    arrange(start,end) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    return()
}
