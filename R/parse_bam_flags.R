#' @title Parse BAM flags
#'
#' @description
#' Read BAM flags as integer values and interpret from binary. Useful for
#' describing a small number of reads in detail. Can return T/F for the presence
#' of flags, but if generating an index for a more complex filter use
#' flag_bam_reads(). If multiple flag values are submitted, a vector of the same
#' length will be returned. If no bam flags are submitted, a table summarizing
#' flag bits is returned.
#'
#' @param bam_flag Integer value to interpret.
#' @param filter_bits Restrict results to specific flags, for instance c(2,3) will only report "Properly aligned" or "Unmapped" flags. If NULL (default), all flags are returned concatenated with a semicolon separator.
#' @param as_boolean Boolean; should the presence of any flag be represented as T/F? Defaults to FALSE.
#'
#' @import dplyr
#' @import magrittr
#' @import tibble
#'
#' @examples
#' parse_bam_flags(return_flag_table=TRUE)
#' parse_bam_flags(1169)
#' parse_bam_flags(1169,filter_bits=11)
#' parse_bam_flags(1169,filter_bits=7,as_boolean=TRUE)
#' parse_bam_flags(1169,filter_bits=8,as_boolean=TRUE)
#'
#' @export

parse_bam_flags <- function(bam_flag,filter_bits=NULL,as_boolean=FALSE){
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter

  tb_flags <- tibble(flag=c(
    "Multiple segments",
    "Properly aligned",
    "Unmapped",
    "Next unmapped",
    "Reverse compliment",
    "Next reverse compliment",
    "First segment",
    "Last segment",
    "Secondary",
    "QC fail",
    "PCR duplicate",
    "Supplementary")) %>%
    mutate(bit = row_number(),
           decimal = 2^c(0:(n()-1)))
  if(missing(bam_flag)){
    res <- tb_flags
  }else{
    if(length(bam_flag) > 1){
      # Apply function to unique values to reduce redundant calls.
      unq_flags <- unique(bam_flag)
      unq_flags <- sapply(unq_flags,parse_bam_flags,filter_bits=filter_bits,as_boolean=as_boolean) %>%
        setNames(as.character(unq_flags))
      res <- unq_flags[as.character(bam_flag)]
    }else{
      if(is.null(filter_bits)) filter_bits <- tb_flags$bit
      res <- tb_flags %>%
        mutate(bit_vals=intToBits(bam_flag)[1:n()],
               report = bit_vals > 0 & bit %in% filter_bits)
      if(as_boolean){
        res <- any(res$report)
      }else{
        res <- res %>%
          filter(report) %>%
          summarize(flags = paste(unlist(flag),collapse=";")) %>%
          unlist(use.names=FALSE)
      }
    }
  }
  return(res)
}
