#' @title
#' Annotate GR
#'
#' @description
#' Given a GRanges object (query) and a second GRanges object to compare it to
#' (subject), identify all overlaps and annotate query by adding a new column
#' (cols_query) with values from a specified column on subject (cols_subject).
#' For instance, if cols_query is set as "peak_overlaps" and cols_subject
#' contains names of each peak in subject, function will return query with a
#' new column called "peak_overlaps" with the associated names from subject. If
#' no overlaps are found for a feature, it is annotated with the value of
#' <na_val>, which defaults to an empty string. If more than one overlap is
#' detected, all values are combined in a list, which can be collapsed into a
#' string separated by semicolons (collapse_as_string; can be slow) or reduced
#' to its first element (first_only). Function can also accept multiple query
#' and subject columns, provided these are the same length.
#'s
#' @param gr_query GRanges object to be annotated.
#' @param gr_subject GRanges object to be overlapped
#' @param cols_query Column(s) to be created in output GRanges object.
#' @param cols_subject Column(s) containing values to be inserted into query columns.
#' @param max_gap Maximum distance between features; defaults to 0bp. Must be set mutually exclusive to min_overlap.
#' @param min_overlap Minimum overlap between features; defaults to 1bp. Must be set mutually exclusive to max_gap.
#' @param as_boolean Should annotations be TRUE/FALSE (overlapped/not overlapped)? Defaults to FALSE. If TRUE, ignores column values to append.
#' @param na_val Value to annotate with in the event no overlaps are detected. Defaults to an empty string.
#' @param collapse_as_string Should annotations from multiple values be collapsed into a single semicolon-separated string? Defaults to FALSE, in which case these are returned as lists.
#' @param first_only Should only the first annotation be returned? Defaults to false; does not occur if <collapse_as_string> is TRUE.
#'
#' @import stringr
#' @import GenomicRanges
#' @import S4Vectors
#' @import tibble
#' @import dplyr
#' @import magrittr
#' @import tidyr
#'
#' @export

annotate_gr     <- function(gr_query,gr_subject,cols_query,cols_subject=NULL,max_gap=-1L,min_overlap=0L,as_boolean=FALSE,na_val = "",collapse_as_string=FALSE,first_only=FALSE){
  olaps   <- findOverlaps(query = gr_query,subject=gr_subject,maxgap = max_gap,minoverlap = min_overlap)

  if(missing(as_boolean)){
    #Default to not boolean.
    as_boolean <- rep(FALSE,length(cols_query))
  }
  if(length(as_boolean) != length(cols_query)){
    #If one as_boolean value provided, recycle it.
    if(length(as_boolean) == 1){
      as_boolean <- rep(as_boolean,length(cols_query))
    }else{
      #If not one, but also different from number of columns, stop.
      stop("If provided, length of 'as_boolean' (",length(as_boolean),
           ") should be either 1 to apply to all columns OR equal to the number of columns (",
           length(cols_query),").")
    }
  }
  for(i in 1:length(cols_query)){
    q_col <- cols_query[i]
    s_col <- cols_subject[i]
    a_bool<- as_boolean[i]

    if(a_bool){
      mcols(gr_query)[q_col] <- FALSE
      mcols(gr_query[queryHits(olaps)])[q_col] <- TRUE
    }else if(is.null(cols_subject)){
      if(is.null(names(gr_subject))){
        names(gr_subject) <- paste0("range_",c(1:length(gr_subject)))
      }
      mcols(gr_query)[q_col] <- na_val
      mcols(gr_query[queryHits(olaps)])[q_col] <- names(gr_subject[subjectHits(olaps)])
    }else{
      if(length(cols_query) != length(cols_subject)){
        stop("Query and subject column names must be the same length.")
      }

      q_tb <- as_tibble(olaps) %>%
        mutate(subject_vals = mcols(gr_subject)[,s_col][subjectHits]) %>%
        group_by(queryHits) %>%
        summarize(subject_vals = list(subject_vals),.groups="drop")
      if(collapse_as_string){
        q_tb <- q_tb %>%
          rowwise %>%
          mutate(subject_vals = paste(subject_vals,collapse=";"))
      }else if(first_only){
        q_tb <- q_tb %>%
          mutate(subject_vals = sapply(subject_vals,function(x) x[[1]]))
      }
      q_tb  <- as.data.frame(q_tb)
      mcols(gr_query)[q_col] <- na_val
      mcols(gr_query)[[q_col]][q_tb$queryHits] <- q_tb$subject_vals
    }
  }
  return(gr_query)
}
