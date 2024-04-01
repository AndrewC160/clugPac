#' @title Trivia stats for basepair distances.
#' 
#' @description
#' Not particualrly useful, but good for trivia/context.
#' 
#' @export

bp_to_distance<- function(bp_value = sum(width(gr_seqs))){
  val_meters  <- bp_value * 0.34 / 1e9
  val_cm      <- val_meters * 100
  
  # Inches.
  val_feet    <- bp_value / 12
  val_miles   <- bp_value / 63360
  val_earths  <- val_miles / 24901
  
  # cm
  val_m       <- bp_value / 100
  val_km      <- val_m / 1000
  
  message(prettyBP(bp_value))
  message("Actual length:")
  message("\t",round(val_meters,digits=3),"m")
  message("\t",round(val_cm,digits=3),"cm")
  message("One inch per bp:")
  message("\t",prettyNumbers(val_feet)," feet")
  message("\t",prettyNumbers(val_miles)," miles")
  message("\t",round(val_earths,digits=2)," equator trips")
}
