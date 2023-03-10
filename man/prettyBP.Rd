% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prettyBP.R
\name{prettyBP}
\alias{prettyBP}
\title{prettyBP}
\usage{
prettyBP(bp_val, force_unit = NULL, digits = 1, sep = "")
}
\arguments{
\item{bp_val}{Numeric value of genomic distance in basepairs}

\item{force_unit}{If provided, unit will be forced. Can be "kb", "bp", "mb", or "MB" for kilobase, basepair, megabase, respectively.}

\item{digits}{Number of digits to be displayed; defaults to 1.}

\item{sep}{Separator between output number and unit, defaults to no space ("100MB").}
}
\description{
Express a number in basepairs with appropriate units.
}
