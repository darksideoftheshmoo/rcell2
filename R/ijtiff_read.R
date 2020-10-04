#' Read TIFF tag information without actually reading the image array.
#'
#' TIFF files contain metadata about images in their _TIFF tags_. This function
#' is for reading this information without reading the actual image.
#'
# @inheritParams read_tif
#' @param path A string. The path to the tiff file to read.
#' @param list_safety A string. This is for type safety of this function. Since
#'   returning a list is unlikely and probably unexpected, the default is to
#'   error. You can instead opt to throw a warning (`list_safety = "warning"`)
#'   or to just return the list quietly (`list_safety = "none"`).
#' @param msg Print an informative message about the image being read?
#' @param frames Which frames do you want to read tags from. Default first frame
#'   only. To read from the 2nd and 7th frames, use `frames = c(2, 7)`, to read
#'   from all frames, use `frames = "all"`.
#'
#' @return A list of lists.
#'
#' @author Simon Urbanek, Kent Johnson, Rory Nolan.
#'
#' @seealso [read_tif()]
#'
#' @examples
#' rcell2::read_tags(path = "data/image_samples/BF_Position001.tif", frames = 1)
#' 
#' @export
read_tags <- function(path, frames = 1, list_safety = "error", msg = TRUE) {
  frames <- prep_frames(frames)
  path <- prep_path(path)
  withr::local_dir(attr(path, "path_dir"))
  if (isTRUE(all.equal(frames, 1,
    check.attributes = FALSE, check.names = FALSE
  ))) {
    return(
      list(frame1 = .Call("read_tags_C", path, 1L, PACKAGE = "rcell2")[[1]])
    )
  }
  tags1 <- read_tags(path, frames = 1)[[1]]
  prep <- prep_read(path, frames, tags1, tags = TRUE)
  out <- .Call("read_tags_C", path, prep$frames, PACKAGE = "rcell2")[prep$back_map]
  frame_nums <- prep$frames[prep$back_map]
  if (!is.na(prep$n_slices) && prep$n_dirs != prep$n_slices) {
    frame_nums <- ceiling(frame_nums / prep$n_ch)
  }
  names(out) <- paste0("frame", filesstrings::nice_nums(frame_nums))
  out
}
