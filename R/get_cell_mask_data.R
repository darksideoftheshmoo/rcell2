#' Read cell masks from BF.out
#'
#' Read boundary and interior pixel data from a BF.out TIFF image.
#' 
#' This function reads cell masks from 16-bit BF.out images onto which boundary
#' and interior pixels have been encoded as intensity values proportional to each
#' object's cellID.
#' 
#' By default the function accounts for the possibility that the image may have
#' a brightfield image background and/or interior pixel intensities offset from
#' boundary pixels. It then first identifies the boundary pixel range, then
#' checks if the image also has offset interior pixels, and finally extracts
#' boundary and interior pixels separately. This behavior can be overridden by
#' explicitly specifying the image type via the \code{blank_bg} and
#' \code{interior_offset} parameters. This may speed up running time.
#'
#' @param path A string. The path to the 16-bit BF.out tiff file to read.
#' @param cell_id_offset the offset with respect to maximum pixel intensity, such
#' that \code{cellID = maximum_intensity - boundary_intensity + cell_id_offset}.
#' @param interior_offset logical. If \code{TRUE} the function expects that
#' cell masks have interior pixels offset from boundary pixels.
#' @param blank_bg logical. If \code{TRUE} the function assumes that the image
#' background (i.e., non-mask pixels) are blank. Set to FALSE if left NULL (the default).
#' 
# @importFrom ijtiff read_tags
#' 
#' @return A data frame of mask pixels ordered by cellID.
#' @export
#'
read_tiff_masks <- function(path, cell_id_offset = -1, interior_offset = NULL, blank_bg = NULL){
  # Set default values for interior_offset and blank_bg, if unspecified by user
  if(is.null(interior_offset)) interior_offset <- TRUE
  if(is.null(blank_bg)) blank_bg <- FALSE
  
  # Read TIFF image
  img <- ijtiff::read_tif(path = path,msg = FALSE)
  
  # Verify 16-bit image and set max image cellID intensity value
  image_bits <- attributes(img)$bits_per_sample
  stopifnot("Expected 16-bit BF.out image." = image_bits==16)
  max_intensity <- 2**image_bits-1+cell_id_offset
  
  # Remove pixels with intensities exceeding max image cellID intensity value.
  # These should be number labels
  img[img>max_intensity] <- 0
  
  # Get sorted vector of unique image intensity values, remove last value (always zero)
  unique_intensity_vals <- sort(unique(as.vector(img)),decreasing=TRUE)
  unique_intensity_vals <- unique_intensity_vals[-length(unique_intensity_vals)]
  
  has_offset <- 0
  
  # If image has BF background and/or if there are interior pixels offset from
  # boundary pixels, we need to find boundary pixel intensity ranges
  if(blank_bg==FALSE || interior_offset==TRUE){
    ## NOTE ##
    # The interior_offset_threshold value is the same as that set in Cell-ID
    # Do not change value here unless also changed in segment.c
    interior_offset_threshold <- 2500
    
    # Get logical index vector indicating which intensity values that are absent from image
    missing_vals <- !(max_intensity:1%in%unique_intensity_vals)
    
    # Loop over missing_vals and count number of consecutive missing intensity
    # values. If this count reaches interior_offset_threshold, then we can be sure
    # that we have looped past all boundary intensity values. The last found
    # intensity value thus belongs to the highest cellID (max_cellid)
    absent_count <- 0
    for(i in 1:length(missing_vals)){
      absent_count <- missing_vals[i]*(absent_count+missing_vals[i])
      if(absent_count==interior_offset_threshold){
        max_cellid <- i-1-interior_offset_threshold
        break
      }
    }
    max_cellid_intensity <- max_intensity - max_cellid
    
    if(interior_offset==TRUE){
      # Calculate offset between boundary and interior pixels for max_cellid objects
      ## NOTE ##
      # This is the same calculation as that done for interior_offset in Cell-ID.
      # Do not change the calculation here unless also changed in segment.c
      interior_intensity_offset <- ((max_cellid+1)%/%interior_offset_threshold+2)*interior_offset_threshold
      
      # Check if offset pixels exist for max_cellid
      if(match(max_cellid_intensity-interior_intensity_offset,unique_intensity_vals,nomatch = 0)){
        max_cellid_offset_intensity <- max_cellid_intensity - interior_intensity_offset
        has_offset <- 1
      }
    }
    
  } else{
    # If background is blank and there are no offset interior pixels, then the
    # smallest intensity value in the image belongs to the highest cellID
    max_cellid_intensity <- unique_intensity_vals[length(unique_intensity_vals)]
  }
  
  # Set threshold for mask pixels
  mask_threshold <- if(has_offset) max_cellid_offset_intensity else max_cellid_intensity
  
  # Get coordinates and values for all mask pixels
  mask_coord <- which(img>=mask_threshold,arr.ind=TRUE)
  mask_vals <- img[mask_coord]
  
  if(has_offset==1){
    # Get boundary indices, coordinates, and values
    boundary_idx <- which(mask_vals>=max_cellid_intensity)
    boundary_coord <- mask_coord[boundary_idx,1:2]
    boundary_vals <- mask_vals[boundary_idx]
    
    # Set boundary intensities to zero after extraction
    mask_vals[mask_vals>=max_cellid_intensity] <- 0
    
    # Get interior indices, coordinates, and values
    interior_idx <- which(mask_vals>=max_cellid_offset_intensity)
    interior_coord <- mask_coord[interior_idx,1:2]
    interior_vals <- mask_vals[interior_idx]
    interior_shift <- max_intensity-max(interior_vals)
    
    # Make mask_data matrix and insert boundary and interior data
    mask_data <- matrix(0,nrow=nrow(mask_coord),ncol=5)
    mask_data[,1] <- c(-boundary_vals,-(interior_vals+interior_shift))+max_intensity
    mask_data[,3:2] <- rbind(boundary_coord,interior_coord)
    mask_data[,4] <- c(boundary_vals,interior_vals)
    
    # Set pixel type (boundary = 1, interior = 2)
    mask_data[,5] <- c(rep(1,nrow(boundary_coord)),rep(2,nrow(interior_coord)))
    
  } else {
    # Make mask_data matrix and insert boundary and interior data
    # Pixel type defaults to 0 (indicating no pixel identify found)
    mask_data <- matrix(0,nrow=nrow(mask_coord),ncol=5)
    mask_data[,1] <- max_intensity-mask_vals
    mask_data[,3:2] <- mask_coord[,1:2]
    mask_data[,4] <- mask_vals
  }
  # Set column names, convert to data.frame and return
  mask_data <- mask_data[order(mask_data[,1]),]
  colnames(mask_data) <- c("cellID","x","y","pix_value","pix_type")
  as.data.frame(mask_data)
}

#' Load cellid boundary data to a dataframe
#' 
#' Using `data_source = "out.tif"` will give **un**ordered boundary pixels.
#' 
#' Using `data_source = "masks.tsv"` will give ordered boundary pixels.
#' 
#' @param data.source Either "out.tif" or "masks.tsv"
#' @param arguments The arguments dataframe used to run cellid (prepared with \code{rcell2::arguments}).
#' @param pixel.type When \code{data.source = "masks.tsv"}, you may choose the pixel "type". At least one of \code{c("i", "b")} for interior and/or boundary pixels ("b" by default).
#' @param close.paths When TRUE and \code{data.source = "masks.tsv"}, append the first row to the end of the data.frame (groping by cellID and pos). Useful for plotting of making closed polygons.
#' @param cell.data When \code{data.source = "out.tif"}, provide the cell.data object (prepared with \code{rcell2::cell.load.alt}).
#' @param tiff.channel When \code{data.source = "out.tif"}, provide the channel name for images holding boundary data ('BF.out' by default).
#' 
#' @return A "cell.boundaries" data.frame.
#' 
#' @export
cell.load.boundaries <- function(data.source,
                                 ...,
                                 arguments = NULL,
                                 pixel.type = "b",
                                 close.paths = FALSE,
                                 cell.data = NULL,
                                 tiff.channel = "BF.out"
                                 ){
  
  option.check <- sum(data.source %in% c("out.tif", "masks.tsv"))
  if(option.check != 1) stop("'data.source' must be either 'out.tif' or 'masks.tsv' ")
  
  if(data.source == "out.tif"){
    
    stopifnot(!is.null(cell.data))
    stopifnot("images" %in% names(cell.data))
    stopifnot(is.data.frame(cell.data$images))
    
    cell.boundaries <-
      cell.data$images %>% 
      filter(channel == tiff.channel) %>% 
      mutate(pos_t.frame = paste(pos, t.frame, sep = "_")) %>% 
      with(setNames(file, pos_t.frame)) %>% 
      lapply(FUN = rcell2::read_tiff_masks, ...) %>% 
      bind_rows(.id = "pos_t.frame") %>% 
      separate(pos_t.frame, into = c("pos", "t.frame"), convert = T)
  }
  
  if(data.source == "masks.tsv"){  
    
    stopifnot(!is.null(arguments))
    stopifnot(is.data.frame(arguments))
    
    cell.boundaries <- 
      rcell2::tsv_paths_from_args(arguments, ...) %>% 
      with(setNames(path, pos)) %>% 
      lapply(readr::read_tsv) %>% 
      bind_rows(.id = "pos") %>% 
      mutate(pos = as.integer(pos),
             # CellID positions are zero-indexed, add 1 and enter the R inferno:
             x=x+1,y=y+1) %>% 
      filter(pixtype %in% pixel.type)
    
    if(close.paths){
      cell.boundaries <- cell.boundaries %>% 
        group_by(pos, cellID) %>% 
        dplyr::slice(c(1:n(), 1))
    }
  }
  
  
  return(cell.boundaries)
}
