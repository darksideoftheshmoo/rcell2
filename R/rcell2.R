#' Yeast Cell Cytometry Suite for CellID in R.
#'
#' It is a rewrite and revamp of the previously awesome Rcell package, by Dr. Alan Bush.
#' Plotting functions from that package have been excluded in the name of minimalism.
#' Thoughfully, \code{\link[ggplot2]{ggplot2}} definitions for all those "cplots" are available in our vignettes (happy copy-pasting!).
#'
#' The full \code{rcell2} functionality is split into three main packages:
#' rcell2, \code{\link[rcell2.cellid]{rcell2.cellid}}, and \code{\link[rcell2.magick]{rcell2.magick}}.
#'
#' @section rcell2 (tidyCell) functions:
#' The tidyCell functions manage CellID's output.
#' It includes functions to load individual cell masks from CellID's identified masks, from generated TSV and/or TIFF files (\link{ijtiff} is required for the latter).
#' Also useful to turn custom data into compatible dataframes for the other functions in this package.
#'
#' @section rcell2.magick functions:
#' \code{\link[rcell2.magick]{magickCell}} renders images from individual cells, based on original images, user defined filters and data from cells.
#' It requires appropriate data frames for observation data (cdata) and for image metadata (paths). It should not be required that the data comes from images processed by CellID. 
#' It relies heavily on the \link{magick} package, so the equirement of the \code{imagemagick} system library might be inconvenient. That is why this awesome feature is optional.
#'
#' \code{\link[rcell2.magick]{shinyCell}} is an R-Shiny based graphical interface, meant to filter, inspect and annotate cells relying on 2D plots of arbitrary variables.
#' It should not be required that the data comes from images processed by CellID, just as \code{\link[rcell2.magick]{magickCell}}.
#' It relies heavily on the \link{magick} package, so the equirement of the \code{imagemagick} system library might be inconvenient. That is why this awesome feature is optional.
#' Also included are two apps: \code{\link[rcell2.magick]{tagCell}} for annotation, and \code{\link[rcell2.magick]{plotApp}} for plotting and subsetting.
#' 
#' @section rcell2.cellid functions:
#' This package includes CellID's source code (written in C). During installation, it compiles a binary executable for Linux, Windows, and Mac systems.
#' The executable is then wrapped seamlessly by the \code{\link[rcell2.cellid]{cell2}} function, such that CellID can now be used programatically within R.
#' If you already have CellID in your system, you are welcome to use it instead of the bundled executable.
#' 
#' @docType package
#' @name rcell2
#' 
#' @seealso \link{ijtiff} \link{magick} \link{shiny}
NULL
