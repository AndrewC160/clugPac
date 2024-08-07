% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/debugDefaultArguments.R
\name{debugDefaultArguments}
\alias{debugDefaultArguments}
\title{Debug default arguments.}
\usage{
debugDefaultArguments(argument_text)
}
\arguments{
\item{argument_text}{String representing all arguments from a function. Can include arguments without a default and/or newlines.}
}
\description{
Quick function for evaluating all default arguments in a given function.
If a function to be tested has at least one default specified in the
function heading, copy the argument and/or the whole set as a string (use
single/double quotes as needed) and pass to debugDefaultArguments().
Function replaces commas with ';' and '=' with '<<-', then uses
parse %>% eval to execute each.

For example, with `arbitrario <- function(input_table,fruit="apple", car="sedan", number=4){...}`, feeding `'input_table,fruit="apple", car="sedan", number=4'` to debugDefaultArguments() will add fruit, car, and number to the top environment with their specified values. Useful for debugging functions with multiple default arguments, but not intended for use in regular code.
}
