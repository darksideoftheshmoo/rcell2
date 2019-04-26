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
