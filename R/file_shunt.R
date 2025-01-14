#' @title File shunt
#'
#' @description
#' Given file names A and B, if file B does not exist copy file A to create B.
#' If file B does exist, ignore file A and return file B. Outputs file B after
#' creating / confirming it exists. Useful for creating local copies of a given
#' file.
#'
#' @param fileA Origin file.
#' @param fileB Destination file.
#' @param create_dir Boolean; should the output directory be created if it is missing? Defaults to FALSE.
#' @param over_write Boolean; should an existing fileB be over-written? Defaults to FALSE.
#'
#' @export

file_shunt <- function(fileA,fileB,create_dir=FALSE,over_write=FALSE){
  if(!file.exists(fileA)) stop("File ",fileA,"not found.")
  if(!file.exists(fileB) | over_write){
    out_dir <- dirname(fileB)
    if(!dir.exists(out_dir)){
      if(!create_dir){
        stop("Output directory not found, create it manually or via 'create_dir=TRUE'.")
      }else{
        dir.create(out_dir,recursive = TRUE)
      }
    }
    file.copy(from = fileA,to = fileB,overwrite = TRUE)
  }
  return(fileB)
}
