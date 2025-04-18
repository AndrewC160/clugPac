% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thin_column.R
\name{thin_column}
\alias{thin_column}
\title{Thin column}
\usage{
thin_column(vals_in)
}
\arguments{
\item{vals_in}{Vector of values to interpret.}
}
\description{
Given a vector of column values, convert any value that is equal to the value
preceding it to empty space. Always returns a character vector, generally
useful for tables to display.

For instance, the column A: "Apples","Apples","Eggs","Eggs","Eggs" will
become "Apples","","Eggs","","".
}
