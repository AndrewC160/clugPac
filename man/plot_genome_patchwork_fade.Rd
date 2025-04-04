% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genome_patchwork_fade.R
\name{plot_genome_patchwork_fade}
\alias{plot_genome_patchwork_fade}
\title{Plot genome patchwork fade.}
\usage{
plot_genome_patchwork_fade(
  gr_regs_in,
  bin_num = 40,
  simplify_labels = TRUE,
  simplify_digits = 2,
  background_fill = "white",
  axis_position = "bottom"
)
}
\arguments{
\item{gr_regs_in}{GRanges object of one or more regions to illustrate.}

\item{bin_num}{Number of bins to split each region into.}

\item{simplify_labels}{Boolean; should genomic values be written in "10MB" annotation? Defaults to false, i.e. "1,000,000".}

\item{simplify_digits}{If <simplify_labels>, how many digits to display when rounding. Defaults to 2.}

\item{background_fill}{Color of the panel background, defaults to "white."}

\item{axis_position}{Position of the x-axis, defaults to "bottom" but can also be "top".}
}
\description{
Given a GRanges object denoting a series of genomic regions on the same
x-axis (for instance, a contig made up of disparate loci), represent it as a
rectangle with fill representing chromosome names. Regions are split into
bins and alpha is used to denote strand information.
}
