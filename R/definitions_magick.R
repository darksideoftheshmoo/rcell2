#' Armar mosaicos cuadrados a partir de un vector de imagenes en magick de cualquier longitud
#' @param images A magick image vector, with images of the same size preferably.
#' @return A single magick image of the squared tile.
# @examples
# square_tile(images)
#' @import magick dplyr
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @export
square_tile <- function(images){
  nRow <- ceiling(sqrt(length(images)))
  nCol <- ceiling(length(images)/nRow)
  
  image.tile <- foreach::foreach(tile_row=0:(nRow-1), .combine=c) %do% {
    
    row_images_index = (1 + tile_row * nCol):(tile_row * nCol + nCol)
    
    # Important for the last row
    if(nRow * nCol > length(images)){
      row_images_index <- row_images_index[row_images_index <= length(images)]
    }
    
    magick::image_append(images[row_images_index])
  }
  
  image.tile <- magick::image_append(image.tile, 
                                     stack = T)
  
  if(length(images) > 50) warning("square_tile says: you have more than 50 images in the input array ¿Are you sure that this is ok?")
  return(image.tile)
}

#' Funcion copada para mostrar fotos basada en magick
#' @details
#' Paths dataframe structure. Output example from \code{glimpse(paths)}:
#' 
#' \preformatted{Columns: 6
#' $ pos     <int> 1
#' $ t.frame <int> 0
#' $ channel <chr> "YFP"
#' $ image   <chr> "picture_filename.tif"
#' $ path    <chr> "/path/to/directory/with/pictures/"
#' $ file    <chr> "/path/to/directory/with/pictures/picture_filename.tif"
#' $ is.out  <lgl> FALSE
#' }
#'
#' @param cdata A Rcell data.frame (not the object).
#' @param paths A paths dataframe with file path, t.frame, position and channel information of each picture.
#' @param max_composite_size Maximum size of the final composite image (this resize is applied last) in pixels. 1000 by default.
#' @param cell_resize Size of the individual cell images. "100x100" by default.
#' @param boxSize Size of the box containing the individual cell images. 50 by default.
#' @param n maximum number of cells to display.
#' @param equalize_images Use magick's function to "equalize" the image when TRUE.
#' @param normalize_images Use magick's function to "normalize" the image when TRUE.
#' @param ch Name of the CellID channel (BF, BF.out, RFP, ...). "BF.out" by default.
#' @param sortVar Variable used to sort the cells. "xpos" by default.
#' @param seed Seed value for sampling of cell images.
#' @param .debug Print more messages if TRUE.
#' @param return_single_imgs If TRUE, return a vector of images instead of a tile.
#' @param  annotation_params a named character vector for magick::annotate options (with one or more of the names "color" "background" "size"). Use NULL to skip annotations.
#' @return A list of two elements: the magick image and the ucids in the image.
# @examples
# magickCell(cdataFiltered, sample_tiff$file, position = sample_tiff$pos, resize_string = "1000x1000")
#' @import magick dplyr
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @export
magickCell <- function(cdata, paths,
                       max_composite_size = 1000, 
                       cell_resize = NULL,
                       boxSize = 50, 
                       n = 100,
                       equalize_images = F, 
                       normalize_images = F,
                       ch = "BF.out",
                       sortVar = "xpos",
                       seed = 1, 
                       .debug=FALSE, 
                       return_single_imgs = FALSE,
                       annotation_params = c(color = "white", 
                                             background = "black")
                       ){
  if(.debug) print("F8")

  # "100x100" pixels
  if(is.null(cell_resize)) cell_resize <- boxSize
  cell_resize_string <- paste0(cell_resize, "x", cell_resize)
  
  # Filter paths by pos and t.frame on cdata
  paths <- filter(paths, pos %in% cdata[,"pos", drop = T] & t.frame %in% cdata[,"t.frame", drop = T])
  
  # Add file path column if absent
  if(!"file" %in% names(paths)){
    paths$file <- normalizePath(paste0(paths$path, "/", paths$image))
    warning("\nWarning (magickCell): 'file' variable absent in paths dataframe. Recreating it from 'path' and 'image' variables.\n")
  }
  
  # Annotation parameters
  updateList <- function(l1, l2){
    common.names <- names(l2)[names(l2) %in% names(l1)]
    l1[common.names] <- l2[common.names]
    return(l1)
  }
  annotation_params_default = c(color = "white",
                                background = "black",
                                size = as.numeric(stringr::str_split(cell_resize_string, "x")[[1]])[1]/7)
  if(!is.null(annotation_params)) annotation_params <- updateList(annotation_params_default, annotation_params)
  

  # Intento con magick
  # magickCell(cdataFiltered, sample_tiff$file, position = sample_tiff$pos, resize_string = "1000x1000")

  getCellGeom <- function(xpos, ypos, boxSize = 50){
    geometry <- magick::geometry_area(width = boxSize,
                                      height = boxSize,
                                      x_off = xpos - base::floor(boxSize/2),
                                      y_off = ypos - base::floor(boxSize/2))
    return(geometry)
  }

  # https://ropensci.org/blog/2017/08/15/magick-10/
  # https://livefreeordichotomize.com/2017/07/18/the-making-of-we-r-ladies/

  n <- min(c(n, nrow(cdata))) # Limit amount of pics to "n"

  # TO-DO
  ## Read only parts of the TIFF matrix to speed up stuff: https://stackoverflow.com/questions/51920316/open-only-part-of-an-image-jpeg-tiff-etc-in-r
  ## Bin the data in 2D and sample 1 cell from within both bin groups.
  ## Can be an animated GIF to show more than 1 cell per bin.
  ## https://rdrr.io/cran/magick/man/image_ggplot.html
  ## funcion cut() https://www.jdatalab.com/data_science_and_data_mining/2017/01/30/data-binning-plot.html
  ## http://finzi.psych.upenn.edu/R/library/ks/html/binning.html
  ## https://rdrr.io/cran/ks/man/binning.html
  ## https://www.rdocumentation.org/packages/npsp/versions/0.7-5/topics/binning

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
    picPath.df <- subset(paths, pos == position & channel %in% ch & t.frame == t_frame)
    picPath <- picPath.df$file[order(match(picPath.df$channel, ch))]  # order paths according to ch argument
    
    stopifnot(length(position) == 1 &length(ucid) == 1 &length(t_frame) == 1) # Checks
    stopifnot(length(picPath) == length(ch) & is.character(picPath)) # Checks

    magick::image_read(picPath) %>%
      {if (equalize_images) magick::image_normalize(.) else .} %>%
      {if (normalize_images) magick::image_equalize(.) else .} %>%
      magick::image_crop(getCellGeom(xpos = cdataSample$xpos[i],
                                     ypos = cdataSample$ypos[i],
                                     boxSize)) %>%
      # Add square black box
      {magick::image_composite(
        magick::image_blank(boxSize, boxSize, "black"),
        ., gravity = "Center"
      )} %>% 
      # Resize
      magick::image_resize(cell_resize_string) %>%
      # Annotate
      {if(is.null(annotation_params)) . else 
        magick::image_annotate(.,
                               text = paste(paste0("Pos", as.character(position)),
                                            paste0("t", t_frame),
                                            ch),
                               color = annotation_params["color"],
                               boxcolor = annotation_params["background"],
                               size = annotation_params["size"],
                               font = "Comic sans",
                               gravity = "SouthEast")} %>%
      {if(is.null(annotation_params)) . else 
        magick::image_annotate(.,
                               text = as.character(ucid),
                               color = annotation_params["color"],
                               boxcolor = annotation_params["background"],
                               size = annotation_params["size"],
                               font = "Comic sans",
                               gravity = "NorthWest")} %>%
      # Add black border
      magick::image_border("black","1x1") %>% 
      # Tile horizontally
      magick::image_append()
    }

  stopifnot(length(imga) == nrow(cdataSample)) # Checks
  
  if(return_single_imgs) return(list("img" = imga,
                                     "ucids" = cdataSample$ucid))

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

  if(nRow*cell_resize >= max_composite_size || nCol*cell_resize >= max_composite_size){
    resize_string <- paste0(max_composite_size, "x", max_composite_size)
    imgb <- magick::image_resize(imgb, resize_string)
  }

  return(list("img" = imgb,
              "ucids" = cdataSample$ucid))
}

#' Funcion copada para mostrar fotos basada en magick
#'
#' @param picPath Image path for the correct position, only one.
#' @param interpolate Passed on to layer params.
#' @return A custom annotation_raster for ggplot.
# @examples
# annotation_magick(.img.path)
#' @import magick grDevices ggplot2
#' @export
annotation_magick <- function(picPath, interpolate = FALSE) {

  # Por ahi era mas facil armarla con NSE: https://wiki.frubox.org/proyectos/atr/rdevel#filter

  stopifnot(length(picPath) == 1 & is.character(picPath)) # Checks

  .img <- picPath %>%
    magick::image_read() %>%
    magick::image_normalize() %>%
    # magick::image_rotate(180) %>%
    magick::image_flip() %>%
    magick::image_fill("none")

  raster <- as.raster(.img)

  .img.info <- magick::image_info(.img)

  raster <- grDevices::as.raster(raster)
  ggplot2::layer(#data = dummy_data(),
    data = data.frame(x = NA),
    mapping = NULL, stat = ggplot2::StatIdentity,
    position = ggplot2::PositionIdentity, geom = ggplot2::GeomRasterAnn, inherit.aes = FALSE,
    params = list(raster = raster,
                  xmin = 0,
                  xmax = .img.info$width,
                  ymin = 0,
                  ymax = .img.info$height,
                  interpolate = interpolate))
}
