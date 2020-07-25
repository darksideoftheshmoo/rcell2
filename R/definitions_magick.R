#' Funcion copada para mostrar fotos basada en magick
#'
#' @param cdata A Rcell data.frame (not the object).
#' @param paths Paths a la imagen de cada posición.
#' @param max_size Maximum size of the final image (this resize is applied last) in pixels. 500 by default.
#' @param cell_resize Size of the individual cell images. "100x100" by default.
#' @param boxSize Size of the box containing the individual cell images. 50 by default.
#' @param n maximum number of cells to display.
#' @param ch Name of the CellID channel (BF, BF.out, RFP, ...). "BF.out" by default.
#' @param sortVar Variable used to sort the cells. "xpos" by default.
#' @param seed Seed value for sampling of cell images.
#' @return Lots of stuff.
# @examples
# magickCell(cdataFiltered, sample_tiff$file, position = sample_tiff$pos, resize_string = "1000x1000")
#' @import magick foreach dplyr
#' @export
magickCell <- function(cdata, paths, 
                       max_size = 500, cell_resize = 100,
                       boxSize = 50, n = 100, 
                       .equalize = F, .normalize = T,
                       ch = "BF.out", sortVar = "xpos", 
                       seed = 1, .debug=FALSE){
  if(.debug) print("F8")

  # "100x100" pixels
  cell_resize_string <- paste0(cell_resize, "x", cell_resize)

  # Intento con magick
  # magickCell(cdataFiltered, sample_tiff$file, position = sample_tiff$pos, resize_string = "1000x1000")

  getCellGeom <- function(xpos, ypos, boxSize = 50){
    geometry <- magick::geometry_area(width = boxSize, 
                                      height = boxSize, 
                                      x_off = xpos - boxSize/2, 
                                      y_off = ypos - boxSize/2)
    return(geometry)
  }

  # https://ropensci.org/blog/2017/08/15/magick-10/
  # https://livefreeordichotomize.com/2017/07/18/the-making-of-we-r-ladies/

  n <- min(c(n, nrow(cdata))) # Limit amount of pics to "n"

  # TO-DO
  # Bin the data in 2D and sample 1 cell from within both bin groups.
  # Can be an animated GIF to show more than 1 cell per bin.
  # https://rdrr.io/cran/magick/man/image_ggplot.html
  # funcion cut() https://www.jdatalab.com/data_science_and_data_mining/2017/01/30/data-binning-plot.html
  # http://finzi.psych.upenn.edu/R/library/ks/html/binning.html
  # https://rdrr.io/cran/ks/man/binning.html
  # https://www.rdocumentation.org/packages/npsp/versions/0.7-5/topics/binning

  # Sample first
  set.seed(seed)
  .cdataSample <- cdata[sample(1:nrow(cdata), n, replace = F),] # sample n rows from cdata
  cdataSample <- .cdataSample[order(.cdataSample[[sortVar]]),  # sort the sample
                       unique(c("pos", "xpos", "ypos", "ucid", "t.frame", sortVar))]  # keep only the necessary columns

  imga <- 
    foreach::foreach(i=1:nrow(cdataSample), .combine=c) %do% {

    position <- cdataSample$pos[i]
    ucid <- cdataSample$ucid[i]
    t_frame <- cdataSample$t.frame[i]
    picPath <- subset(paths, pos == position & channel == ch & t.frame == t_frame)$file
    
    stopifnot(length(position) == 1 &length(ucid) == 1 &length(t_frame) == 1) # Checks
    stopifnot(length(picPath) == 1 & is.character(picPath)) # Checks
    
    magick::image_read(picPath) %>%
      {if (.normalize) magick::image_normalize(.) else .} %>% 
      {if (.equalize) magick::image_equalize(.) else .} %>% 
      magick::image_crop(getCellGeom(xpos = cdataSample$xpos[i],
                                     ypos = cdataSample$ypos[i],
                                     boxSize)) %>%
      magick::image_resize(cell_resize_string) %>%
      magick::image_annotate(text = paste(paste0("Pos", as.character(position)), 
                                          paste0("t", t_frame),
                                          ch),
                             size = as.numeric(stringr::str_split(cell_resize_string, "x")[[1]])[1]/7,
                             color = "white",
                             boxcolor = "black",
                             font = "Comic sans",
                             gravity = "SouthEast") %>%
      magick::image_annotate(text = as.character(ucid),
                             # text = paste0(as.character(ucid), "t", t_frame),
                             size = as.numeric(stringr::str_split(cell_resize_string, "x")[[1]])[1]/7,
                             color = "white",
                             boxcolor = "black",
                             font = "Comic sans",
                             gravity = "NorthWest") %>%
      magick::image_border("black","1x1")
    }
  
  stopifnot(length(imga) == nrow(cdataSample)) # Checks

  nRow <- ceiling(sqrt(n))
  nCol <- ceiling(n/nRow)

  imgb <- foreach::foreach(i=0:(nRow-1), .combine=c) %do% {

    j = (1 + i*nCol):(i*nCol + nCol)
    j = j[j <= n]

    magick::image_apply(imga[j], function(i){
      magick::image_blank(boxSize, boxSize, "black") %>%
        magick::image_resize(cell_resize_string) %>%
        magick::image_composite(i) }) %>%
      magick::image_append()
  }
  imgb <- magick::image_append(imgb, stack = T)

  if(nRow*cell_resize >= max_size || nCol*cell_resize >= max_size){
    resize_string <- paste0(max_size, "x", max_size)
    imgb <- magick::image_resize(imgb, resize_string)
  }

  return(list("img" = imgb,
              "ucids" = cdataSample$ucid))
}

