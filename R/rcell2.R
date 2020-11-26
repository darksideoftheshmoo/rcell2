#' Yeast Cell Cytometry Suite for CellID in R.
#'
#' It is a rewrite and revamp of the previously awesome Rcell package, by Dr. Alan Bush.
#' Plotting functions from that package have been excluded in the name of minimalism.
#' Thoughfully, \link{ggplot2} definitions for all those "cplots" are available in our vignettes (happy copy-pasting!).
#'
#' The cellMagick package provides three categories of important functions:
#' cellMagick, tidyCell and shinyCell.
#'
#' @section tidyCell functions:
#' The tidyCell functions run CellID and/or manage it's output.
#' Also useful to turn custom data into compatible dataframes for the other functions in this package.
#'
#' @section cellMagick functions:
#' Renders images from individual cells, based on original images, user defined filters and data from cells.
#' It should not be required that the data comes from images processed by CellID.
#' Requirement of the imagemagick library might bother some, so this awesome feature is optional.
#'
#' @section shinyCell functions:
#' R-Shiny based graphical interface to filter cells by arbitrary variables, and inspect and annotate cells.
#' It should not be required that the data comes from images processed by CellID.
#'
#' @docType package
#' @name rcell2
#' 
#' @seealso \link{magick}
NULL
