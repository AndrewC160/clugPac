#' @title Flag BAM reads
#'
#' @description
#' Interpret BAM flag integers and check for the presence of specific flags.
#' Returns a Boolean vector the same length as the input flag vector. A custom
#' logic function can be applied, as well: logic_function = all will require all
#' flags to be found. Defaults to any, in which case at least one flag must be
#' found. In the event multiple flags are submitted, unique flag values will be
#' interpreted separately and returned as a corresponding vector of the same
#' length as the input.
#'
#' Bit flags:
#'                            bit decimal
#' Multiple segments           1       1
#' Properly aligned            2       2
#' Unmapped                    3       4
#' Next unmapped               4       8
#' Reverse compliment          5      16
#' Next reverse compliment     6      32
#' First segment               7      64
#' Last segment                8     128
#' Secondary                   9     256
#' QC fail                    10     512
#' PCR duplicate              11    1024
#' Supplementary              12    2048
#'
#' @param bam_flags Integer value(s) to interpret.
#' @param filter_bits Restrict results to specific flags, for instance c(2,3) will only report "Properly aligned" or "Unmapped" flags.
#' @param logic_function Function for returning T/F when parsing flags; defaults to any.
#'
#' @examples
#' #Check for reverse-compliment reads.
#' flag_bam_reads(1169,5)
#' flag_bam_reads(129,5)
#'
#' #Identify reverse-compliment reads that are ALSO the last segment.
#' flag_bam_reads(bam_flags=c(369,353),c(5,9),logic_function=any)
#'
#' @export

flag_bam_reads <- function(bam_flags,filter_bits=11,logic_function=any){
  if(length(bam_flags)>1){
    # Apply function to unique flag values to avoid redundant calls.
    flag_vals <- unique(bam_flags)
    unq_vals  <- sapply(flag_vals,flag_bam_reads,filter_bits=filter_bits,logic_function=logic_function) %>%
      setNames(as.character(flag_vals))
    unq_vals[as.character(bam_flags)]
  }else{
    logic_function(filter_bits %in% which(intToBits(bam_flags) > 0))
  }
}
