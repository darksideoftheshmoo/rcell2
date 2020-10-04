#' @useDynLib rcell2
#' @importFrom magrittr '%>%' '%<>%' '%T>%'
#' @importFrom rlang '%||%'
#' @importFrom zeallot '%<-%'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "n_slices", "n_ch", "ij_n_ch", "n_imgs"))
}

.onUnload <- function(libpath) {
  library.dynam.unload("rcell2", libpath)
}
