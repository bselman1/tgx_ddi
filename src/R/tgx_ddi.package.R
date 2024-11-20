#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end

#' Define a package wide environment for storing package variables
this_pkg <- new.env(parent = emptyenv())

#' @title
#' Helper function to load a file defined by this package in the *inst* directory.
get_pkg_filepath = function(...){
  system.file(..., package = "tgxddi", mustWork = TRUE)
}

#' @title Check if an object is an error
#' @description Check to see if the provided object is an S3 class that inherits from "try-error"
#' @param x An object to check
is.error = function(x) inherits(x, "try-error")

get_error_msg = function(e) {
  condition = attr(e, "condition")
  if (is.null(condition)) {
    return("")
  }
  return(conditionMessage(condition))
}

NULL