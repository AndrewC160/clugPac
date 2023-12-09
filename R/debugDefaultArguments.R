#' @title
#' Debug default arguments.
#'
#' @description
#' Quick function for evaluating all default arguments in a given function.
#' If a function to be tested has at least one default specified in the
#' function heading, copy the argument and/or the whole set as a string (use
#' single/double quotes as needed) and pass to debugDefaultArguments().
#' Function replaces commas with ';' and '=' with '<<-', then uses
#' parse %>% eval to execute each.
#'
#' For example, with `arbitrario <- function(input_table,fruit="apple", car="sedan", number=4){...}`, feeding `'input_table,fruit="apple", car="sedan", number=4'` to debugDefaultArguments() will add fruit, car, and number to the top environment with their specified values. Useful for debugging functions with multiple default arguments, but not intended for use in regular code.
#'
#' @param argument_text String representing all arguments from a function. Can include arguments without a default and/or newlines.
#'
#' @import stringr
#'
#' @export

debugDefaultArguments <- function(argument_text){
  argument_text %>%
    gsub(",",";",.) %>%
    gsub("=","<<-",.) %>%
    parse(text=.) %>%
    eval
}
#argument_text <- c("a,z,x,y=4,g=2")

#debugDefaultArguments2(argument_text)


# debugDefaultArguments <- function(argument_text){
#   arg_list  <- gsub("\n","",argument_text)
#   #From https://stackoverflow.com/questions/35347537/using-strsplit-in-r-ignoring-anything-in-parentheses
#   arg_list  <- strsplit(argument_text, '\\([^)]+,(*SKIP)(*FAIL)|,\\s*', perl=TRUE)[[1]]
#   arg_list  <- str_match(arg_list,pattern = "(.+)=(.+)")
#   arg_list  <- arg_list[complete.cases(arg_list),]
#   apply(arg_list,MARGIN = 1, function(x) eval(parse(text=paste(x[2],"<<-",x[3]))))
# }
