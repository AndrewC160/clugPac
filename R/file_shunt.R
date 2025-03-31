#' @title File shunt
#'
#' @description
#' Given file names A and B, if file B does not exist copy file A to create B.
#' If file B does exist, ignore file A and return file B. Outputs file B after
#' creating / confirming it exists. Useful for creating local copies of a given
#' file. Can provide a conditional BASH command (i.e., returns empty text if
#' fileB exists) as well, which is useful for producing basic scripts that can
#' be executed outside of R.
#'
#' @param fileA Origin file.
#' @param fileB Destination file.
#' @param create_dir Boolean; should the output directory be created if it is missing? Defaults to FALSE.
#' @param over_write Boolean; should an existing fileB be over-written? Defaults to FALSE.
#' @param dry_run Boolean; should file not actually be copied? Defaults to FALSE, otherwise will print a message and return the filename.
#' @param output_bash_cmd Boolean; should a BASH command be returned that would accomplish the file copy? Defaults to FALSE, otherwise the text output is a copy command that can be run in BASH.
#'
#' @export

file_shunt <- function(fileA,fileB,create_dir=TRUE,over_write=FALSE,dry_run=FALSE,output_bash_cmd=FALSE){
  if(!file.exists(fileA)) stop("File ",fileA,"not found.")
  out_cmd   <- ""
  if(!file.exists(fileB) | over_write){
    out_dir <- dirname(fileB)
    if(!dir.exists(out_dir)){
      if(!create_dir){
        stop("Output directory ",dirname(fileB)," not found, create it manually or via 'create_dir=TRUE'.")
      }else{
        out_cmd <- paste("cp -p",fileA,fileB)
        dir.create(out_dir,recursive = TRUE)
      }
    }else{
      out_cmd   <- paste("cp",fileA,fileB)
    }
    if(dry_run){
      message("Will copy ",fileA," to ",fileB,".")
    }else if(!output_bash_cmd){
      file.copy(from = fileA,to = fileB,overwrite = TRUE)
    }
  }
  if(output_bash_cmd){
    out_txt <- out_cmd
  }else{
    out_txt <- fileB
  }
  return(out_txt)
}
