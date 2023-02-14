#' @title
#' bigBedToBed
#'
#' @description
#' Wrapper for UCSC tools; convert bigbed format file to BED format given an
#' input genomic range of interest.
#'
#' @param gr_in GRange to query.
#' @param bb_file Filename of bigBed format file.
#' @param exe_location Location of UCSC tools' bigBedToBed executable.
#'
#' @import tibble
#' @import dplyr
#' @import data.table
#' @import GenomicRanges
#'

#DOES NOT WORK YET,DO NOT EXPORT.

gr_in   <- GRanges("chr8",IRanges(102000000,102350000))
bigBedToBed   <- function(gr_in,bb_file,exe_location="/asclab/projects/shared/modules/ucsc-tools/ucsc-tools-20170321/bin/bigBedToBed",max_ranges=10){
  if(missing(gr_in)) stop("No ranges of interest provided; exiting.")
  if(missing(bb_file)) stop("No bigBed file provided; exiting.")
  if(!file.exists(bb_file)) stop("bigBed file '",bb_file,"' not found; exiting.")
  if(length(gr_in) > 1){
    if(length(gr_in) > max_ranges){
      stop(paste0("Do you mean to query ",length(gr_in)," ranges from ",basename(bb_file),"? If so increase max_ranges."))
    }
    lapply(1:length(gr_in), function(i){
      bigBedToBed(gr_in = gr_in[i],
                  bb_file = bb_file,
                  exe_location = exe_location)
    })
  }
  #bigBedToBed -chrom=chr2 -start=0 -end=100000 ./k24.Unique.Mappability.bb stdout
  seq_nm    <- as.character(seqnames(gr_in))
  start_val <- start(gr_in)
  end_val   <- end(gr_in)
  tb_map  <-
    paste0(exe_location,
         " -chrom=",seq_nm,
         " -start=",start_val,
         " -end=",end_val,
         " ",bb_file,
         " stdout") %>%
    fread(cmd=.,sep = "\t")


}
