#' Extract info from an ImageJ-style `TIFFTAG_DESCRIPTION`.
#'
#' @inheritParams prep_read
#'
#' @return A named list with elements `n_imgs`, `n_slices`, `ij_n_ch`, `n_ch`.
#'
#' @noRd
translate_ij_description <- function(tags1) {
  n_imgs <- NA_integer_
  n_slices <- NA_integer_
  ij_n_ch <- FALSE
  n_ch <- tags1$samples_per_pixel %||% 1
  if ("description" %in% names(tags1) &&
      startsWith(tags1$description, "ImageJ")) {
    ij_description <- tags1$description
    if (stringr::str_detect(ij_description, "channels=")) {
      n_ch <- filesstrings::first_number_after_first(
        ij_description,
        "channels="
      )
      ij_n_ch <- TRUE
    }
    n_imgs <- filesstrings::first_number_after_first(ij_description, "images=")
    n_slices <- calculate_n_slices(ij_description)
    if ((!is.na(n_slices) && !is.na(n_imgs)) &&
        ij_n_ch &&
        n_imgs != n_ch * n_slices) {
      stop(
        "
        The ImageJ-written image you're trying to read says in its
        TIFFTAG_DESCRIPTION that it has {n_imgs} images of
        {n_slices} slices of {n_ch} channels. However, with {n_slices}
        slices of {n_ch} channels, one would expect there to be
        {n_slices} x {n_ch} = {n_ch * n_slices} images.
        ", "
        This discrepancy means that the `ijtiff` package can't read your
        image correctly.
        ", "
        One possible source of this kind of error is that your image
        is temporal and volumetric. `ijtiff` can handle either
        time-based or volumetric stacks, but not both."
      )
    }
  }
  list(n_imgs = n_imgs, n_slices = n_slices, ij_n_ch = ij_n_ch, n_ch = n_ch)
}


#' Check that the `frames` argument has been passed correctly.
#'
#' @param frames An integerish vector. The requested frames.
#'
#' @return `TRUE` invisibly if everything is OK. The function errors otherwise.
#'
#' @noRd
prep_frames <- function(frames) {
  checkmate::assert(
    checkmate::check_string(frames),
    checkmate::check_integerish(frames, lower = 1)
  )
  if (is.character(frames)) {
    frames %<>% tolower()
    if (!startsWith("all", frames)) {
      custom_stop(
        "If `frames` is a string, it must be 'all'.",
        "You have `frames = '{frames}'`."
      )
    }
    frames <- "all"
  }
  frames
}

#' Get information necessary for reading the image.
#'
#' While doing so, perform a check to see if the requested frames exist.
#'
#' @param path The path to the TIFF file.
#' @param frames An integerish vector. The requested frames.
#' @param tags1 The tags from the first image (directory) in the TIFF file. The
#'   first element of an output from [read_tags()].
#' @param tags Are we prepping the read of just tags (`TRUE`) or an image
#'   (`FALSE`).
#'
#' @return A list with seven elements.
#' * `frames` is the adjusted frame numbers (allowing for _ImageJ_  stuff),
#'  unique and sorted.
#' * `back_map` is a mapping from `frames` back to its non-unique, unsorted
#'  original; that would be `frames[back_map]`.
#'  * `n_ch` is the number of channels.
#'  * `n_dirs` is the number of directories in the TIFF image.
#'  * `n_slices` is the number of slices in the TIFF file. For most, this is the
#'   same as `n_dirs` but for ImageJ-written images it can be different.
#'  * `n_imgs` is the number of images according to the ImageJ
#'   `TIFFTAG_DESCRIPTION`. If not specified, it's `NA`.
#'  * `ij_n_ch` is `TRUE` if the number of channels was specified in the ImageJ
#'   `TIFFTAG_DESCRIPTION`, otherwise `FALSE`.
#'
#' @noRd
prep_read <- function(path, frames, tags1, tags = FALSE) {
  frames <- prep_framesframes
  frames_max <- max(frames)
  
  # c(n_imgs, n_slices, ij_n_ch, n_ch) %<-% translate_ij_description(
  #   tags1
  # )[c("n_imgs", "n_slices", "ij_n_ch", "n_ch")]
  translated_ij_description <- translate_ij_description(tags1)
  n_imgs <- translated_ij_description["n_imgs"]
  n_slices <- translated_ij_description["n_slices"]
  ij_n_ch <- translated_ij_description["ij_n_ch"]
  n_ch <- translated_ij_description["n_ch"]
  
  path <- prep_path(path)
  withr::local_dir(attr(path, "path_dir"))
  n_dirs <- .Call("count_directories_C", path, PACKAGE = "ijtiff")
  if (!is.na(n_slices)) {
    if (frames[[1]] == "all") {
      frames <- seq_len(n_slices)
      frames_max <- n_slices
    }
    if (frames_max > n_slices) {
      custom_stop("
      You have requested frame number {frames_max} but there
      are only {n_slices} frames in total.
                ")
    }
    if (ij_n_ch) {
      if (n_dirs != n_slices) {
        if (!is.na(n_imgs) && n_dirs != n_imgs) {
          custom_stop(
            "
          If TIFFTAG_DESCRIPTION specifies the number of images, this must be
          equal to the number of directories in the TIFF file.
          ",
            "Your TIFF file has {n_dirs} directories.",
            "Its TIFFTAG_DESCRIPTION indicates that it holds {n_imgs} images."
          )
        }
        if (tags) {
          frames <- frames * n_ch - (n_ch - 1)
        } else {
          frames <- purrr::map(
            frames * n_ch,
            ~ .x - rev((seq_len(n_ch) - 1))
          ) %>%
            unlist()
        }
      }
    }
  } else {
    if (frames[[1]] == "all") {
      frames <- seq_len(n_dirs)
      frames_max <- n_dirs
    }
    if (frames_max > n_dirs) {
      custom_stop("
      You have requested frame number {frames_max} but there
      are only {n_dirs} frames in total.
                ")
    }
  }
  good_frames <- sort(unique(frames))
  back_map <- match(frames, good_frames)
  list(
    frames = as.integer(good_frames),
    back_map = back_map,
    n_ch = n_ch,
    n_dirs = n_dirs,
    n_slices = ifelse(is.na(n_slices), n_dirs, n_slices),
    n_imgs = n_imgs,
    ij_n_ch = ij_n_ch
  )
}

#' Prepare the path to a TIFF file for a function that will read from that file.
#'
#' The [fs::path_file()] is returned. The calling function is expected to call
#' `withr::local_dir(fs::path_dir())`.
#'
#' @param path A string. The path to a TIFF file.
#'
#' @return A string. The [fs::path_file()]. This has an attribute `path_dir`
#'   with the path to be passed to [withr::local_dir()].
#'
#' @noRd
prep_path <- function(path) {
  checkmate::assert_string(path)
  path %<>% stringr::str_replace_all(stringr::coll("\\"), "/") # windows safe
  checkmate::assert_file_exists(path)
  structure(fs::path_file(path), path_dir = fs::path_dir(path))
}
