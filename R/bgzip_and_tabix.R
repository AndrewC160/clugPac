#' @title
#' BGzip and Tabix
#'
#' @description
#' Accepts an input bedGraph or BED file for bgzipping and tabix indexing.
#' Executable locations for Tabix, Bedtools, and BGzip are hard-coded.
#'
#' @param file_in Input file, must not be gzipped. If a .gz ending is present, the file will be presumed to not yet be compressed.
#' @param silent Boolean, should messages denoting steps be suppressed? Defaults to FALSE.
#' @param bedtools_exe BedTools executable.
#' @param bgzip_exe BGzip executable.
#' @param tabix_exe Tabix executable.
#'
#' @import tictoc
#'
#' @export

bgzip_and_tabix   <- function(file_in,silent=FALSE,
                              bedtools_exe="~/src/bedtools2/bin/bedtools",
                              bgzip_exe="/asclab/projects/shared/modules/samtools/samtools-1.2/bin/el7/bgzip",
                              tabix_exe="/asclab/projects/shared/modules/samtools/samtools-1.2/bin/el7/tabix"){
  if(grepl(".gz$",file_in)){
    if(!silent) message("Removing '.gz' ending...")
    new_file_in   <- gsub(".gz$","",file_in)
    file.rename(file_in,new_file_in)
    file_in <- new_file_in
  }
  if(!silent) tic()
  if(!silent) message("Sorting...")
  system(paste0("grep -v 'track' ",file_in," |",bedtools_exe," sort -i - > ",file_in,".srtd"))
  if(!silent) message("BGzipping...")
  system(paste0(bgzip_exe," -c ",file_in,".srtd > ",file_in,".gz"))
  if(!silent) message("Building Tabix index...")
  system(paste0(tabix_exe," -p bed ",file_in,".gz"))
  file.remove(file_in)
  file.remove(paste0(file_in,".srtd"))
  if(!silent) toc()
  return(paste0(file_in,".gz"))
}
