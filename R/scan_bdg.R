#' @title
#' Scan Tabix-indexed BedGraph file.
#'
#' @description
#' Scan one or more BGzipped and Tabix-indexed files within regions specified
#' in <gr_regions> GRange object. Function is run recursively if more than one
#' bgzipped file is specified, and results are concatenated into a single
#' table. Uses Rsamtools::TabixFile for scanning and returns a tibble
#'
#' @param bgz_files List of filenames of tabix-indexed BDG files.
#' @param name_vec Vector of names for 'name' column. If missing, names can be retrieved from bgz_files, otherwise will use file base name. Must be the same length as bgz_files.
#' @param gr_regions GenomicRanges containing regions of interest.
#' @param col_names Column names of BedGraph file. Defaults to standard 5 columns (seq,start,end,score,name), but can be changed to include additional columns if necessary.
#'
#' @import Rsamtools
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @import GenomicRanges
#'
#' @export

scan_bdg        <- function(bgz_files,name_vec,gr_regions,col_names){
  if(missing(col_names)){
    col_names   <- c("seqnames","start","end","score","region")
  }
  if(length(bgz_files) > 1){
    if(is.null(names(bgz_files)) & missing(name_vec)){
      names(bgz_files)  <- paste0("bdg_",c(1:length(bgz_files)))
    }else{
      if(length(name_vec) != length(bgz_files)) stop("Differing number of files/names provided.")
      names(bgz_files)  <- name_vec
    }
    tb_out  <- lapply(names(bgz_files),function(nm) {
      scan_bdg(bgz_files = bgz_files[[nm]],gr_region = gr_regions,col_names=col_names) %>%
        as_tibble %>%
        mutate(name = nm)
      }) %>%
      do.call(rbind,.)
  }else{
    tb_file <- Rsamtools::TabixFile(file = bgz_files)
    txt_lst <- Rsamtools::scanTabix(tb_file,param=gr_regions)
    if(length(txt_lst) == 0){
      return(NULL)
    }
    if(missing(name_vec)){
      name_value  <- gsub(".bdg.gz","",basename(tb_file$path))
    }else{
      name_value  <- name_vec
    }
    tb_out  <- lapply(names(txt_lst), function(nm) {
      txt   <- txt_lst[[nm]]
      paste0(txt,"\t",nm)
    }) %>%
      do.call(c,.)
    if(length(tb_out) == 1){
      #Add newline in cases where only one value returned.
      tb_out <- paste0(tb_out,"\n")
    }
    tb_out  <- fread(text=tb_out,sep="\t",col.names = col_names) %>%
      as_tibble %>%
      mutate(name = name_value)
  }
  return(tb_out)
}
