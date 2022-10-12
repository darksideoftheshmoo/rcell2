# rcell2

Functions to analyze [Cell-ID](https://github.com/darksideoftheshmoo/cellID-linux)'s single-cell cytometry data in R, using a tidy framework.

rcell2's functionality is split into four packages:

* The main rcell2 package offers functions to load Cell-ID's output to data.frames, and image manipulation based on EBImage. A development version of this package is available in the [`rcell.dev`](https://github.com/darksideoftheshmoo/rcell2/tree/rcell2.dev) branch.
* Cell-ID, the image segmentation software, has been wrapped in the [`rcell.cellid`](https://github.com/darksideoftheshmoo/rcell2/tree/rcell2.cellid) package. It offers functions to run CellID from R, and an rmarkdown template showcasing advanced functionality.
* The cell tiling and graphic filtering apps, built on R-Shiny and [magick](https://github.com/ropensci/magick), are available in the [`rcell.magick`](https://github.com/darksideoftheshmoo/rcell2/tree/rcell2.magick) package.
* The [`rcell2.examples`](https://github.com/darksideoftheshmoo/rcell2.examples) package contains notebooks on general usage, and on several classification and analysis methods.

This package is very well tested in baker's yeast data, and R version 4+.

## Preview

Provide a _defocused_ brightfield image to CellID, and _voila_:

<img src="doc/TFP.png" width="350"> <img src="doc/TFP.out.png" width="350">

> Segmentation of yeast cells in a single position.

The image is segmented, cells are identified and tracked over time, and features are computed from morphology and fluorescent signal distribution.

<img src="doc/preview_tile.png" width="700">

> Time series images of one cell, showing different acquisition channels. 

With Rcell2, you can load an analize the CellID results freely, using standard R packages.

<img src="doc/preview_trace.png" width="700">

> Background corrected fluorescent signal concentration VS time, plotted with ggplot2.

# Installation

## R Dependencies

Most of the dependencies are listed in the `DESCRIPTION` file, and should install automatically.

We suggest installing the [```tidyverse```][1] meta-package (and use it, you'll not regret it) and [```devtools```][2]:

```r
install.packages(c("tidyverse", "devtools"))
```

In addition, install [```EBImage```][3] package (required to look at cells) by copying and running the following script:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EBImage")
```


<!-- DEPRECATED SECTION: magick was moved to its own package.

Also install `imagemagick` on your system; this is required by R's `magick` package. All the major operating systems are supported by ImageMagick: https://imagemagick.org/script/download.php

For Ubuntu, Arch and macOS (use [homebrew](https://brew.sh/)!!!)  these commands may come in handy:

```sh
# Aptitude
sudo apt install imagemagick

# Pacman
sudo pacman -S imagemagick

```

For macOS see: https://imagemagick.org/script/download.php

## CellID dependency (optional)

The `cell2` function can be used to run CellID directly from R.

It supports newer CellID functions, available at the repo's `mask_mod` branch: https://github.com/darksideoftheshmoo/cellID-linux/tree/mask_mod

Older CellID versions may work, as long as the new options are not enabled.

Visit the CellID repository to get installation instructions.

-->


## Installing the package

Install using `remotes`, directly from github repo:

```
# rcell2 package (dev branch)
remotes::install_github("darksideoftheshmoo/rcell2", ref = "rcell2.dev")

# cellid package
remotes::install_github("darksideoftheshmoo/rcell2-cellid")

# shiny-magick package
remotes::install_github("darksideoftheshmoo/rcell2-magick")
```

Example data and analysis notebooks can be found in the examples package: [`rcell2.examples`](https://github.com/darksideoftheshmoo/rcell2.examples/tree/main).

# New tools

## CellID bundled in R

> This feature is available in the new "rcell2.cellid" package.

* This package bundles, compiles and wraps our improved Cell-ID binary.
* Save, load and manipulate cell masks and boundaries within R.
* Run it in parallel automagically, backed by R's `foreach`.
* Automate scanning parameter ranges, by running Cell-ID programatically.

## R-Shiny tools for cytometry data

> This feature is available in the new "shiny-magick" package.

<!-- DEPRECATED SECTION: magick was moved to its own package.

An R-Shiny app will help users filter data graphically, with live image previews.
While the image part may be tailored for data from fluorescence microscopy experiments, the graphical filter in this app is general purpose (i.e. useful in standard cell cytometry).

There is also another small app to "tag" single cells in the dataset with user defined options.

New functions based on the `magick` package help build image tiles and strips of individual cells.

See:

* `?rcell2::magickCell()`
* `?rcell2::cellStrips()`
* `?rcell2::cellSpreadPlot()`
* `?rcell2::shinyCell()`
* `?rcell2::tagCell()`
* `?rcell2::plotApp()`

-->

Single-cell visualization and image plotting tools:

* Reads single cells images of all channels into `magick` objects.
* Use any of `magick`'s functions to manipulate your images.
* Creates tiles and mosaics of single cells, time courses, of one or many imaging/fluorescence channels.
* Plot of a 2D grid of "representative" single cell images in a scatterplot (similar to EBImage). Implemented in functions: cellSpread, cellSpreadPlot, and "Pics" type plot in shinyCell.


## Hu Moment functions for raw cell segmentation data

We implemented the `Hu moments` descriptors in R, and use them on masks generated by CellID. Note that the masks must be generated by the the CellID [`mask_mod` branch](https://github.com/darksideoftheshmoo/cellID-linux/tree/mask_mod) either by TSV output or by encoding CellIDs in the pixel intensities of boundary and/or interior points.

We recommend using the [`rcell.cellid`](https://github.com/darksideoftheshmoo/rcell2/tree/rcell2.cellid) package to generate the required input.


## K-means filtering functions

The kmeans algotrithm helps filter cells based on clustering of CellID's variables computed from morphological and fluorescence information.

Use k-means and check out images of cells in each cluster. Then, filter them easily by cluster number.

See: `?rcell2::kmeans_clustering`


## Development notebooks

Example image datasets and analysis notebooks can be found in the rcell2 examples companion package: [`rcell2.examples`](https://github.com/darksideoftheshmoo/rcell2.examples/tree/main).

The `extdata/testings` directory holds many rmarkdown notebooks, where we explore different analysis approaches to single cell images and cytometry:

* Spatial distribution of fluorescent signals.
* Pattern detection in time-series.
* Classification examples.
* Cell boundary curvature analysis and alignment.
* ...


# Todo

* ~~Plot of a 2D grid of "representative" single cell images in a scatterplot (similar to EBImage).~~ Implemented in rcell2.magick`: cellSpread, cellSpreadPlot, and "Pics" type plot in shinyCell.
* ~~Generic shiny function for filtering points in a custom ggplot~~. Implemented in rcell2.magick`: plotApp.
* 1D / histogram filtering support in the filtering app.
* Per-facet filtering in shinyCell.
* cellStrips support in shinyCell image viewer.

!['Automating' comes from the roots 'auto-' meaning 'self-', and 'mating', meaning 'screwing'.](https://imgs.xkcd.com/comics/automation.png)

[1]:https://www.tidyverse.org/
[2]:https://github.com/r-lib/devtools
[3]:https://bioconductor.org/packages/release/bioc/html/EBImage.html
