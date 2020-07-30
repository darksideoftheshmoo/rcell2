# rcell2
functions to analyze Cell-ID data in a tidy framework


# Installation

## dependencies

Install the [```tidyverse```][1] meta-package (and use it, you'll not regret it) and [```devtools```][2]. In addition, install [```EBImage```][3] package (required to look at cells) by copying and running the following script. 



```
install.packages("tidyverse", "devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EBImage")
```

## installing the package

Install using devtools, from github repo:

```
devtools::install_github("gerbeldo/rcell2")

```

# rcellid & cellMagick

This package bundles the CellID image segmentation program, and has utility functions to generate single-cell data from microscopy images.

After importing data with Rcell2, an R-Shiny app will help users filter data graphically, with live image previews.

While the image part may be tailored for data from fluorescence microscopy experiments, the graphical filter in this app is general purpose (i.e. useful in standard cell cytometry).

## Todo

* Plot of a 2D grid of "representative" single cell images in a scatterplot.
* Generic shiny function for filtering points in a custom ggplot.
* 1D / histogram filtering support.

## Installation

1. Download this pacakge.
2. Install external library dependencies `libtiff` and `glib-2.0`.
3. Install the package from source.

### Installing dependencies

For Ubuntu, Arch and Brew (macOS):

```
# Aptitude
sudo apt install libglib2.0-dev libtiff5-dev

# Pacman
sudo pacman -S glib2 libtiff

# Brew
brew install glib
brew install libtiff
```

### Installing the pakage from source

```
install.packages("devtools")

path_to_file <- "cellMagick_0.1.3.tar.gz"

devtools::install_local(path_to_file, dependencies = TRUE)
```

## Future work

This package will include the CellID program source, and wrap it in a single R function. For now it can only call a CellID executable provided by you, dear user.

[1]:https://www.tidyverse.org/
[2]:https://github.com/r-lib/devtools
[3]:https://bioconductor.org/packages/release/bioc/html/EBImage.html
