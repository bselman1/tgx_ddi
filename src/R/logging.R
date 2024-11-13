#' @title Get package logfile path.
#' @description
#' Get the value of the LOG_FILEPATH variable stored in the package environment
#' if available. Will return NULL if the variable has not yet been set.
get_logfile_path = function() {
  if (!exists("LOG_FILEPATH", envir = this_pkg)) {
    return(NULL)
  } else {
    return(get("LOG_FILEPATH", envir = this_pkg))
  }
}


#' @title Set package logfile path.
#' @description
#' Set the value of the LOG_FILEPATH variable stored in the package environment.
#' If this value is set, all logged messages will be sent to the file specified 
#' as well as the console.
#' @export
set_logfile_path = function(logfile_path = NULL) {
  assign("LOG_FILEPATH", value = logfile_path, envir = this_pkg)
}


#' @title Log a message
#' @description
#' Helper method that will log a message to the logfile (if configured) and send
#' echo the message to the console using the [base::message] function
log_msg <- function(msg, file) {
  message(msg);
  
  logfile_path = get_logfile_path()
  if (!is.null(logfile_path)) {
    cat(msg, "\n", file = logfile_path, append = TRUE);
  }
}