
#' An S4 class to represent hubviz model fitting results.
#'
#' @slot data data
#' @slot init initial value
#' @slot result result

setClass( Class="hubviz",
    representation=representation(
        data="matrix",
        init="list",
        result="list"
        )
)

