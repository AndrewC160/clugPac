#' @title Combine PNGs
#'
#' @description
#' Combine a list of PNG image filenames into a single PNG raster. Note that aspect ratios of input
#' files must be respected (no GGPlot2 rescaling takes place once a PNG is read as a raster). For
#' this reason, to combine two plots side-by-side in a 10" x 5" output plot, both input files should
#' first be saved as 5" x 5" PNG files.
#'
#' @param p_file_list List of PNG files to combine.
#' @param n_row Number of rows for the output image. Defaults to 1.
#' @param n_col Number of columns for the output image. If not provided, <n_col> is calculated based on the number of plots provided and <n_row>.
#'
#' @import png
#'
#' @export

combine_pngs <- function(p_file_list,n_row=1,n_col){
  if(missing(n_col)) n_col <- ceiling(length(p_file_list)/n_row)
  px <- lapply(p_file_list,function(png_fl){
    as.raster(readPNG(png_fl)) %>% rasterGrob(interpolate=FALSE)
  })
  marrangeGrob(grobs=px,nrow=n_row,ncol=n_col,top=NULL)
}
