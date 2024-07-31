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
#' 
#' @export

file_shunt <- function(fileA,fileB){
  if(!file.exists(fileB)){
    file.copy(fileA,fileB)
  }
  return(fileB)
}
