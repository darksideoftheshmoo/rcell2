% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/definitions_magick.R
\name{square_tile}
\alias{square_tile}
\title{Armar mosaicos cuadrados a partir de un vector de imagenes en magick de cualquier longitud}
\usage{
square_tile(images, nRow = NULL, nCol = NULL)
}
\arguments{
\item{images}{A magick image vector, with images of the same size preferably.}

\item{nRow}{An integer indicating number of rows. If specified, nCol must be specified too.}

\item{nCol}{An integer indicating number of rows. If specified, nRow must be specified too.}
}
\value{
A single magick image of the squared tile.
}
\description{
Armar mosaicos cuadrados a partir de un vector de imagenes en magick de cualquier longitud
}
