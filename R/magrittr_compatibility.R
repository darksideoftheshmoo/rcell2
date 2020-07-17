# #' filter.cell.data
# #'
# #' @param .data cell.data object
# #' @param ... extra arguments for dplyr
# #'
# #' @return filtered tibble ()
# #' @export
# #'
# #' @examples
# filter.cell.data <- function(.data, ...) {
#     l <- rlang::enquos(...)
#     dplyr::filter(.data$data, !!! l)
# }
#


#' d
#'
#' convenience function to use when piping.
#'
#' d(X) is equivalent to X$data
#'
#' @param cell.data cell.data object
#'
#' @return a tibble
#' @export
#'
#' @examples
d <- function(cell.data) {
    cell.data[['data']]
}
