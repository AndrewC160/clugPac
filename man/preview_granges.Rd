% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preview_granges.R
\name{preview_granges}
\alias{preview_granges}
\title{Preview GRanges.}
\usage{
preview_granges(list_grs, window_coords = NULL)
}
\arguments{
\item{list_grs}{List of GenomicRanges, either named or not.}

\item{window_coords}{GenomicRanges object representing window to focus on. Defaults to NULL.}
}
\description{
Given a list of GenomicRange objects, this function produces a summary plot
using GGplot2, focused on whatever genomic region was provided in
"window_coords". If no window_coords provided, attempts to find a region
with the most overlaps between all input ranges.
}
