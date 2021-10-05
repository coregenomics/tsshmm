## usethis namespace: start
#' @useDynLib tsshmm, .registration = TRUE
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
    library.dynam.unload("tsshmm", libpath) # nocov
}
