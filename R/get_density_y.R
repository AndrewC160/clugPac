#' @title Get density Y-values.
#'
#' @description Given a column of specific values, determine the approximate y-
#' position of each value in a kernel density plot (i.e. geom_density() in
#' ggplot2). Useful for adding labels to density plots.
#'
#' @param vals_in Values to approximate.
#' @param density_vals All values used to measure the kernel density.
#'
#' @import stats
#' @import tibble
#' @import dplyr
#' @import magrittr
#'
#' @export

get_density_y <- function(vals_in,density_vals){
  mutate <- dplyr::mutate
  filter <- dplyr::filter
  select <- dplyr::select

  ob_dens <- density(density_vals)
  tb  <- tibble(x=ob_dens$x,y=ob_dens$y)
  get_closest_y <- function(x_val){
    tb %>%
      mutate(dif = abs(x - x_val)) %>%
      filter(dif == min(dif)) %>%
      filter(row_number() == 1) %>% #In the off chance that the min. value occurs twice.
      select(y) %>%
      unlist %>%
      return()
  }
  sapply(vals_in,get_closest_y) %>%
    setNames(NULL) %>%
    return()
}
