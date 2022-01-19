#' Run ImageJ FFT filter macro from R
#' 
#' Passes an ImageJ FFT filter on files matching \code{BF.*.tif} in the target directory.
#' 
#' The modified images are saved to a new subdirectory with default name: "filtered/".
#' 
#' @param pic.path Path to the directory containing the image files.
#' @param imagej.path Path to the ImageJ binary (a path to "ImageJ-linux64" or equivalent).
#' @param script.path Path to the ImageJ macro. Defaults to built-in macro.
#' @export
#' 
imagej.fft.filter <- function(
  pic.path,
  imagej.path = "~/Software/Fiji.app/ImageJ-linux64",
  script.path = system.file("inst/FFT_filter_on_BFs_R.txt", package = "rcell2")){
  
  is.dir <- file.info(pic.path)$isdir
  
  if(is.dir) {
    pic.path <- paste0(normalizePath(pic.path), "/")
  } else {
    stop("The provided 'pic.path' is not a directory.")
  }
  
  command <- paste(
    normalizePath(imagej.path),
    "--headless -macro",
    normalizePath(script.path),
    pic.path
  )
  
  system(command)
}
