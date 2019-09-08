#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("vdjbaseVis version: ",packageVersion(pkgname)))
  invisible()
}
