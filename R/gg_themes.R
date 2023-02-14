#' @title
#' Clug themes
#'
#' @description
#' Assorted GGplot2 themes, can be used similarly to 'theme_minimal()'.
#' Additional arguments to theme() can also be supplied.
#'
#' @param theme_name Name of theme of interest ("genome","scatter","barplot").
#'
#' @import ggplot2
#'
#' @export

gg_theme <- function(theme_name = "genome",...){
  theme_name  <- tolower(theme_name)
  if(theme_name == "genome"){
    t_out   <- theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(size=0.25,color="lightgray",linetype = "dotted"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank()
    )
  }else if(theme_name == "scatter"){
    t_out   <- theme(
      panel.background = element_blank(),
      panel.border = element_rect(size=0.5,fill=NA,color="black"),
      panel.grid.major = element_line(size=0.25,color="lightgray"),
      strip.background = element_blank()
    )
  }else if(theme_name == "barplot"){
    t_out   <- theme(
      panel.background = element_blank(),
      panel.border = element_rect(size=0.5,fill=NA,color="black"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size=0.25,color="lightgray"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle=45,hjust=1,vjust=1),
      strip.background = element_blank()
    )
  }else{
    stop("No theme named '",theme_name,"' found...")
  }
  t_out   <- t_out + theme(...)
  return(t_out)
}
