#' Run ImageJ FFT filter macro from R
#' 
#' Passes an ImageJ FFT filter on files matching \code{BF.*.tif} in the target directory.
#' 
#' The modified images are saved to a new subdirectory with default name: "filtered/".
#' 
#' @param pic.path Path to the directory containing the image files.
#' @inheritParams imagej.macro.run
#' @export
#' 
imagej.fft.filter <- function(
  pic.path,
  script.path = system.file("inst/imagej_macros/FFT_filter_on_BFs_R.txt", 
                            package = "rcell2")){
  
  is.dir <- file.info(pic.path)$isdir
  
  if(is.dir) {
    pic.path <- paste0(normalizePath(pic.path), "/")
  } else {
    stop("The provided 'pic.path' is not a directory.")
  }
  
  imagej.macro.run(script.path, pic.path)
}

# #' Run ImageJ open on a path
# #' 
# #' @param pic.path Path to the image file.
# #' @inheritParams imagej.macro.run.headless
# #' @export
# #' 
# imagej.open <- function(
#   pic.path,
#   script.path = system.file("inst/imagej_macros/open.ijm", 
#                             package = "rcell2"),
#   wait = F){
#   
#   is.dir <- file.info(pic.path)$isdir
#   
#   if(!is.dir) {
#     pic.path <- paste0(normalizePath(pic.path))
#   } else {
#     stop("The provided 'pic.path' is a directory.")
#   }
#   
#   imagej.macro.run(script.path, pic.path)
# }


#' Run headless ImageJ Macro file
#' 
#' @param imagej.path Path to the ImageJ binary (a path to "ImageJ-linux64" or equivalent).
#' @param script.path Path to the ImageJ macro. Defaults to built-in macro.
#' @inheritParams base::system
#' 
imagej.macro.run <- function(
  script.path,
  ...,
  imagej.path = "~/Software/Fiji.app/ImageJ-linux64",
  wait = T, headless = T){
  command <- paste(
    normalizePath(imagej.path),
    {if (headless) "--headless -macro" else "-macro"},
    normalizePath(script.path),
    ...
  )
  
  base::system(command, wait = wait)
}

