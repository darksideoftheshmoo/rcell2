#' calculate number of time frames a cell appears in
#' 
#' In case a cell was filtered out, its n.tot will be 0.
#'
#' @param data data.frame
#'
#' @return data.frame with n.tot variable
#' @export
#'
calculate_ntot <- function(data) {
  
  counts <- data %>% 
    dplyr::select(-dplyr::contains("n.tot"), ucid, qc) %>% 
    dplyr::filter(qc) %>% 
    dplyr::count(ucid)
  
  data %>% 
    dplyr::left_join(counts, by = "ucid") %>% 
    tidyr::replace_na(list(n = 0)) %>% 
    dplyr::rename(n.tot = n)
  
}


#' Update number of frames a cell appears in
#'
#' @param object cell.data
#'
#' @return cell.data
#' @export
#'
update_ntot <- function(object) {
  
  stopifnot("This is not a timecourse experiment" = length(unique(object$data$t.frame)) > 1)
  
  object$data <- calculate_ntot(object$data)
  object$variables$all <- union(object$variables$all, "n.tot")
  
  return(object)
}


