#' @title Isometrify facets
#'
#' @description
#' Given a table with very basic categorical faceting and x/y values, return a
#' a table columns for creating a pseudo-3D "isometric" view which
#' emulates three dimensions by shifting facet tracks in the x and y directions
#' and ordering geoms such geoms in the "foreground" will obstruct geoms in the
#' "background." See also: Diablo games.
#'
#' Function is very brute-force; it returns columns that can be attached to an
#' original table using the "row_num" index column, x and y values, their
#' shifted minimums (x/y_shift), and their shifted values (x/y_shifted). The
#' depth column is generated using the factor levels of the <values_facet>
#' argument; if it is not a factor it will be made into one using the values in
#' the order submitted.
#'
#' @param values_facet Vector of values containing facet information (i.e. categories to be spread in the Z-direction).
#' @param values_x Vector of x values.
#' @param values_y Vector of y values.
#' @param x_rng_in Vector of length 2 defining the x-range limits. If not provided, defaults to the range of <values_x>.
#' @param y_rng_in Vector of length 2 defining the y-range limits. If not provided, defaults to the range of <values_y>.
#' @param x_rng_expand Fraction of x-axis to expand when simulating the shift. All z-facet x-shifts will sum up to this value. Defaults to 0.25.
#' @param y_rng_expand Fraction of y-axis to expand when simulating the shift. All z-facet y-shifts will sum up to this value. Defaults to 0.25.
#' @param include_demo Should a demonstration isometric plot be provided? Defaults to FALSE.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @import ggplot2
#'
#' @export

isometrify_table <- function(values_facet,values_x,values_y,
                             x_rng_in=NULL,y_rng_in=NULL,
                             x_rng_expand=0.25,y_rng_expand=0.25,
                             include_demo=FALSE){
  if(is.null(x_rng_in)) x_rng <- range(values_x)
  if(is.null(y_rng_in)) y_rng <- range(values_y)
  x_rng_shift_max <- x_rng_expand * diff(x_rng)
  y_rng_shift_max <- y_rng_expand * diff(y_rng)
  x_rng_shift <- x_rng_shift_max / length(unique(values_facet))
  y_rng_shift <- y_rng_shift_max / length(unique(values_facet))

  tb <- tibble(facet=values_facet,
               x = values_x,
               y = values_y,
               row_num = c(1:length(values_facet)))

  if(!is.factor(tb$facet)) tb <- mutate(tb,facet=factor(facet,levels=select(.,facet) %>% unlist %>% unique))

  facet_levs <- levels(tb$facet)

  shift_x <- function(x_in,depth_in=length(levels(tb$facet))-1){
    x_in + depth_in * x_rng_shift
  }
  shift_y <- function(y_in,depth_in=length(levels(tb$facet))-1){
    y_in + depth_in * y_rng_shift
  }

  tb_p <- tb %>%
    mutate(depth = as.integer(facet)-1,
           x_shift = depth * x_rng_shift,
           y_shift = depth * y_rng_shift,
           x_shifted = x_shift + x,
           y_shifted = y_shift + y,
           facet = factor(facet,levels=rev(facet_levs)))

  tb_base <- tibble(idx=c(1:4),
                    x=c(x_rng[1],
                        x_rng[2],
                        x_rng[2] %>% shift_x,
                        x_rng[1] %>% shift_x),
                    y=c(y_rng[1],
                        y_rng[1],
                        y_rng[1] %>% shift_y,
                        y_rng[1] %>% shift_y),
                    geom="base")

  p <- ggplot(tb_p,aes(x=x,y=y,group=facet)) + geom_line()
  pb <- ggplot_build(p)

  x_brks <- pb$layout$panel_scales_x[[1]]$break_positions() %>% na.omit %>% as.double
  y_brks <- pb$layout$panel_scales_y[[1]]$break_positions() %>% na.omit %>% as.double

  tb_grid <- rbind(
    tibble(x=unique(c(x_brks,x_rng[2])),
           y=y_rng[1]) %>%
      mutate(xend=shift_x(x),
             yend=shift_y(y),
             set="bottom"),
    tibble(x=c(x_rng[1],x_brks,x_rng[2]) %>% unique %>% shift_x,
           y=shift_y(y_rng[1])) %>%
      mutate(xend = x,
             yend=shift_y(y_rng[2]),
             set="x_back"),
    tibble(x=x_rng[1],y=y_rng[1],
           xend=x_rng[1],
           yend=y_rng[2],
           set="x_left"),
    tibble(x=x_rng[1],
           y=y_rng[2],
           xend=shift_x(x_rng[1]),
           yend=shift_y(y_rng[2]),
           set="top_left"),
    tibble(x=shift_x(x_rng[1]),
           y=shift_y(y_rng[2]),
           xend=shift_x(x_rng[2]),
           yend=shift_y(y_rng[2]),
           set="top_back"),
    tibble(x=x_rng[1],
           y=c(y_rng[1],y_brks,y_rng[2]) %>% unique,
           xend=shift_x(x_rng[1])) %>%
      mutate(yend=shift_y(y),
             set="y_left"),
    tibble(x=shift_x(x_rng[1]),
           y=c(y_rng[1],y_brks,y_rng[2]) %>% unique %>% shift_y,
           xend=shift_x(x_rng[2]),
           yend=y,
           set="y_back"),
    tibble(x=x_rng[1],
           y=y_rng[1],
           xend=x_rng[2],
           yend=y_rng[1],
           set="x_front")
  )
  list_out <-
    list(data=tb_p,
         grid=tb_grid,
         x_brks=x_brks,
         y_brks=y_brks,
         shift_y_func=shift_y,
         shift_x_func=shift_x)

  if(include_demo){
    list_out$example_plot <-
      ggplot(tb_p,aes(x=x_shifted,ymin=y_shift,ymax=y_shifted,fill=facet)) +
        scale_x_continuous(expand=c(0,0),breaks=x_brks) +
        scale_y_continuous(expand=c(0,0),breaks=y_brks) +
        geom_segment(tb_grid,mapping=aes(x=x,y=y,xend=xend,yend=yend),
                     color="grey30",
                     inherit.aes=FALSE) +
        geom_ribbon() +
        theme_minimal() +
        theme(panel.grid = element_blank())
  }
  return(list_out)
}
