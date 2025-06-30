#' @title SessionInfoSorted
#' 
#' @description
#' Returns sessionInfo() output after sorting base, other, and non-attached package lists.
#' 
#' @export

sessionInfoSorted   <- function(){
  si <- sessionInfo()
  order1 <- order(names(si$otherPkgs))
  order2 <- order(names(si$loadedOnly))
  si$basePkgs   <- sort(si$basePkgs)
  si$otherPkgs  <- si$otherPkgs[order1]
  si$loadedOnly <- si$loadedOnly[order2]
  return(si)
}