#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onLoad <- function(libname, pkgname) {
  packageStartupMessage(paste0("vdjbaseVis version: ",packageVersion(pkgname)))
  invisible()
}
