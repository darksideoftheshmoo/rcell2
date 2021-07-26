#' Run ImageJ FFT filter macro from R
#' 
#' Replaces files matching \code{BF.*.tif} in the target directory.
#' 
#' @param pic.path Path to the directory containing the image files.
#' @param imagej.path Path to the ImageJ binary (a path to "ImageJ-linux64" or equivalent).
#' @param script.path Path to the ImageJ macro. Defaults to built-in macro.
#' 
imagej.fft.filter <- function(
  pic.path,
  imagej.path = "~/Software/Fiji.app/ImageJ-linux64",
  # script.path = "~/Projects/Academia/Doctorado/gitlabs_acl/scripts_confocal_vic/FFT_filter_on_BFs_R.txt"){
  script.path = system.file("inst/FFT_filter_on_BFs_R.txt", package = "rcell2")){
  system(paste(
    normalizePath(imagej.path),
    "--headless -macro",
    normalizePath(script.path),
    paste0(dirname(pic.path), "/")
  ))
}
