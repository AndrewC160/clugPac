#' @title
#' prettyTitle
#'
#' @description
#' Given a string, function capitolizes either the first or every word (separated
#' by 'sep'). Runs recursively if given a list or vector, and if the input is a
#' factor with levels and 'preserve_levels' is TRUE, existing levels are also
#' changed and then reapplied. Otherwise, vectors are returned as character type.
#'
#' @param argument_text String representing all arguments from a function. Can include arguments without a default and/or newlines.
#'
#' @import stringr
#' @import stringi
#'
#' @export

prettyTitle <- function(str_in,sep="_",caps="first",preserve_levels=TRUE,preserve_ackronyms=TRUE){
  #Functions.
  cap_func <- function(txt_in = "somatic"){
    paste0(toupper(substr(txt_in,1,1)),substr(txt_in,2,nchar(txt_in)))
  }
  camel_func <- function(txt_in = "CamelCase",save_acks=TRUE){
    #Parse CamelCase and snakeCase, and send to lowercase, but respect ackronyms.
    c_txt   <- unlist(str_split(txt_in,pattern=""))

    #First letter is always boundary capitol (i.e, even in snakeCase).
    cap_pos <- c(1,grep("[A-Z]",unlist(str_split(txt_in,pattern = "")),ignore.case = FALSE)) %>% unique

    #Ignore subsequent capitols ("DNABaseLevel" becomes "DNA base level")
    is_first<- c(TRUE,!unlist(sapply(seq_along(cap_pos), function(i) cap_pos[i] == cap_pos[i-1] + 1)))
    is_last <- c(na.omit(unlist(sapply(seq_along(cap_pos), function(i) cap_pos[i+1] > cap_pos[i] + 1))),TRUE)
    is_term <- cap_pos == nchar(txt_in)
    #Boundary: Capitol OR first capitol in ackronym.
    bounds  <- cap_pos[is_first | (is_last & !is_term)]

    #Do not lowercase ackronyms (ackronym: previous letter was capitol end of
    #string, UNLESS last capitol is not the last letter in the string).
    is_ack  <- c(na.omit(sapply(seq_along(cap_pos), function(i) !is_first[i+1])),FALSE) | is_term
    acks    <- cap_pos[is_ack]
    if(any(acks) & save_acks){
      c_txt[-acks] <- tolower(c_txt[-acks])
    }else{
      c_txt <- tolower(c_txt)
    }

    sapply(seq_along(bounds), function(i) {
      w_strt <- as.integer(bounds[i])
      w_end  <- as.integer(bounds[i + 1]) - 1
      if(is.na(w_end)){
        w_end<- nchar(txt_in)
      }
      return(paste(c_txt[w_strt:w_end],collapse=""))
    })
  }

  if(length(str_in) > 1){
    txt <- sapply(str_in,prettyTitle,sep=sep,caps=caps,preserve_levels=preserve_levels)
    if(is.factor(str_in) & preserve_levels){
      lv_nms  <- sapply(levels(str_in),prettyTitle,sep=sep,caps=caps)
      txt <- factor(txt,levels=lv_nms)
    }
    return(txt)
  }
  if(sep=="camel"){
    txt <- camel_func(str_in,save_acks = preserve_ackronyms)
  }else{
    txt <- unlist(str_split(str_in,pattern=sep,simplify = FALSE))
  }
  if(caps == "first"){
    txt <- paste(cap_func(txt[1]),paste(txt[-1],collapse=" "),sep=" ")
  }else if(caps == "all"){
    txt <- paste(sapply(txt,cap_func),collapse=" ")
  }else{
    stop("Argument 'caps' can be either 'first' or 'all'.")
  }
  return(trimws(txt))
}
