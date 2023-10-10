# Original functions ####

#' Read cell masks from BF.out
#'
#' Read boundary and interior pixel data from a BF.out TIFF image.
#' 
#' This function reads cell masks from 16-bit "BF.out" images produced by Cell-ID,
#' onto which boundary and interior pixels have been encoded as intensity values
#'  proportional to each segmented object's cellID.
#' 
#' By default the function accounts for the possibility that the image may have
#' a brightfield image background and/or interior pixel intensities offset from
#' boundary pixels. It then first identifies the boundary pixel range, then
#' checks if the image also has offset interior pixels, and finally extracts
#' boundary and interior pixels separately. This behavior can be overridden by
#' explicitly specifying the image type via the \code{blank_bg} and
#' \code{interior_offset} parameters. This may speed up running time.
#' 
#' Limitations: this function only uses information from 16-bit "brightfield" 
#' images output by CellID. The fluoresence channels output images are 8-bit,
#' and do not contain the necessary information. This might be relevant if the
#' fluorescent images were aligned to the BF by CellID, such that the pixel
#' coordinates differ between brightfield and fluorescence output images.
#' If this is a problem for you, then use \code{cell.load.boundaries} with
#' \code{data.source = 'masks.tsv'} to read CellID's TSV output files (for 
#' that see \code{cell2}).
#' 
#' See CellID's mask_mod branch for more details: https://github.com/darksideoftheshmoo/cellID-linux/blob/mask_mod/README.md#branch-notes
#'
#' @param path A string. The path to the 16-bit BF.out tiff file to read.
#' @param cell_id_offset integer. The pixel intensity offset with respect to maximum pixel intensity, such
#' that \code{cellID = maximum_intensity - boundary_intensity + cell_id_offset}.
#' @param interior_offset logical. If \code{TRUE} the function expects that
#' cell masks have interior pixels offset from boundary pixels (set to TRUE if NULL; the default).
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
  
  # Check dependency  
  if(!requireNamespace("ijtiff")){
    stop("read_tiff_masks requires the 'ijtiff' package, which is not installed by default. Please install it to use this function :)")
  }
  
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


#' Generate data.frame with paths to coordinate TSV files
#' 
#' Convenience function to generate input for \code{hues_from_tsv_files2}.
#'
#' @param arguments An "arguments" dataframe, as produced by \code{rcell2.cellid::arguments}.
#' @param position_pattern A regex pattern with one group for the integer position number, extracted from the directory name holding the TSV file.
#' @param tsv_pattern A regex pattern matching the name of the TSV file (\code{readr::read_tsv} supports reading "gz" compressed files direclty).
#' @export
tsv_paths_from_args <- function(arguments,
                                position_pattern = ".*Position(\\d+)$",
                                tsv_pattern = "^out_all_masks.tsv(?:\\.gz)?$"){
  
  tsv_paths <- dir(unique(arguments$output), 
                   pattern = tsv_pattern, 
                   full.names = T)
  
  if(length(tsv_paths) != length(unique(arguments$output))) 
    stop("tsv_paths_from_args error: the number of TSVs found does not match the arguments length")
  
  tsv_files_df <- 
    data.frame(
      pos = as.integer(sub(position_pattern, "\\1", dirname(tsv_paths))),
      path = tsv_paths
    )
  
  return(tsv_files_df)
}

#' Load CellID boundary data to a dataframe
#' 
#' @details 
#' 
#' Using \code{data_source = "out.tif"} will give **un** ordered boundary pixels.
#' See \code{read_tiff_masks} for details and caveats.
#' CellID must have been run with the "encode_cellID_in_pixels" option enabled.
#' 
#' Using \code{data_source = "masks.tsv"} will give ordered boundary pixels.
#' CellID must have been run with the "output_coords_to_tsv" option enabled.
#' However, in some cases, the pixel sequence may not be "manifold"
#' (see: https://github.com/darksideoftheshmoo/cellID-linux/issues/11).
#' Also, do not expect them to be consistent in their direction (clockwise or 
#' anti-clockwise), and double check if the pixels are aligned to the output images
#' (see: https://github.com/darksideoftheshmoo/rcell2/issues/29).
#' 
#' @param data.source Either "out.tif" or "masks.tsv".
#' @param arguments The arguments dataframe used to run cellid (prepared with \code{rcell2.cellid::arguments}).
#' @param pixel.type When \code{data.source = "masks.tsv"}, you may choose the pixel "type". At least one of \code{c("i", "b")} for interior and/or boundary pixels ("b" by default).
#' @param flags Used to subset the input files by CellID's "flag" field. Each flag corresponds to an imaging channel, according to a "mapping" found in \code{cell.data$mapping} (available when using get_cell_data).
#' @param close.paths When TRUE and \code{data.source = "masks.tsv"}, append the first row to the end of the data.frame (groping by cellID and pos). Useful for plotting of making closed polygons.
#' @param cell.data When \code{data.source = "out.tif"}, you **must** provide the cell.data object (prepared with \code{rcell2.cellid::get_cell_data}).
#' @param tiff.channel When \code{data.source = "out.tif"}, provide the channel name for images holding boundary data ('BF.out' by default).
#' @inheritParams tsv_paths_from_args
#' @inheritParams read_tiff_masks
#' @return A "cell.boundaries" data.frame.
#' 
#' @export
cell.load.boundaries <- function(data.source,
                                 # tsv_paths_from_args args
                                 position_pattern = ".*Position(\\d+)$",
                                 tsv_pattern = "^out_all_masks.tsv(?:\\.gz)?$",
                                 # read_tiff_masks args
                                 cell_id_offset = -1, interior_offset = NULL, blank_bg = NULL,
                                 # Other args
                                 arguments = NULL,
                                 pixel.type = "b",
                                 flags = NULL,
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
      lapply(FUN = read_tiff_masks, 
             cell_id_offset = cell_id_offset, 
             interior_offset = interior_offset, 
             blank_bg = blank_bg) %>% 
      bind_rows(.id = "pos_t.frame") %>% 
      separate(pos_t.frame, into = c("pos", "t.frame"), convert = T)
  }
  
  if(data.source == "masks.tsv"){
    
    stopifnot(!is.null(arguments))
    stopifnot(is.data.frame(arguments))
    
    cell.boundaries <- 
      tsv_paths_from_args(arguments, 
                          position_pattern=position_pattern,
                          tsv_pattern=tsv_pattern) %>% 
      with(setNames(path, pos)) %>% 
      lapply(readr::read_tsv,
             # col_types = "iiiiif") %>%  # types: 5 int columns, 1 factor column 
             col_types = readr::cols(
               cellID=readr::col_integer(),
               t.frame=readr::col_integer(),
               flag=readr::col_integer(),
               x=readr::col_integer(),
               y=readr::col_integer(),
               pixtype=readr::col_factor(levels = c("b", "i"))
             )) %>%  # types: 5 int columns, 1 factor column 
      bind_rows(.id = "pos") %>% 
      mutate(pos = as.integer(pos),
             # CellID positions are zero-indexed, add 1 and enter the R inferno:
             x=x+1,y=y+1) %>% 
      filter(pixtype %in% pixel.type)
    
    if(!is.null(flags)){
      cell.boundaries <- cell.boundaries %>% filter(flag == flags)
    }
    
    if(close.paths){
      cell.boundaries <- cell.boundaries %>% 
        group_by(pos, cellID) %>% 
        dplyr::slice(c(1:n(), 1))
    }
  }
  
  
  return(cell.boundaries)
}

# TIFF masks section ####

#' Load cell masks from BF.out TIFF files
#'
#' It is imperative that the TIFF images contain only cell boundaries, with pixel intensity values such that \code{CellID = (2^image_bits -1) - pixel_value}.
#' 
#' Such TIFF images are automatically generated by the mask_mod branch of the CellID program found at: https://github.com/darksideoftheshmoo/cellID-linux/tree/mask_mod
#' 
#' @param images An images dataframe, from  \code{rcell2::load_cell_data} or \code{rcell2.cellid::get_cell_data}.
#' @param image_bits an integer indicating the bit-depth on the TIFF images, such that maximum intensity is \code{image_bits^2 -1} and minimum is zero. Typically 8 or 16.
#' @param cell_id_offset the offset respect to maximum pixel intensity, such that \code{cellID = maximum_intensity - boundary_intensity + cell_id_offset}.
#' @param return_points if TRUE it will add a "masks" dataframe to the cell_data object, containing the mask coordinates.
#' @export
#' 
load_tiff_masks <- function(images, image_bits, cell_id_offset = -1, return_points = F){
  paths <- images %>% 
    mutate(file = paste0(path, "/", image)) %>% 
    filter(channel == "BF.out")
  
  tiff.masks.df <- apply(paths, MARGIN = 1, FUN = function(pic_metadata){
    # Extract xy coordinates list from the BF.out mask tiff and derive cellID from the intensities
    pic_df <- mask_df_from_tiff(tiff_path = pic_metadata["file"],
                                image_bits = image_bits, 
                                cell_id_offset = cell_id_offset)
    
    return(pic_df)
  })
  
  masks_df <- bind_rows(tiff.masks.df)
  return(masks_df)
  
}

#' Read cell masks from BF.out
#' 
#' Generate a dataframe with boundary coordinates from mask tiff file.
#' 
#' This function is alternative to Andy's much smarter \code{read_tiff_masks} function.
#' 
# @importFrom tiff readTIFF
#' @keywords internal
#' 
mask_df_from_tiff <- function(tiff_path, image_bits, cell_id_offset = -1){
  
  # Check dependency
  if(!requireNamespace("tiff")){
    stop("mask_df_from_tiff: requires functions from the 'tiff' package, which is not installed by default.")
  }
  
  # Read masks tiff
  pic <- tiff::readTIFF(tiff_path)
  
  # Convert intensity value to 16-bit integers
  pic <- pic * (2^image_bits - 1)
  
  # Replace zeros with NA
  pic[pic == 0] <- NA
  
  # Convert pic matrix to dataframe of coordinates per mask/value
  pic_df <- data.frame(y=rep(1:nrow(pic), ncol(pic)), 
                       x=rep(1:ncol(pic), each = nrow(pic)),
                       pix_value=c(pic))
  # Clear NA rows
  pic_df <- pic_df[!is.na(pic_df[["pix_value"]]),]
  
  # Convert integer intensity value to CellID
  pic_df$cellID <- factor((2^image_bits -1) - pic_df[["pix_value"]] + cell_id_offset)
  
  return(pic_df)
}

# TSV masks section ####

#' Generate a data.frame with XY mask data from CellID's TSV files
#' 
#' @details Do not specify masks_tsv_path, it is generated by tsv_paths_from_dir automatically.
#' 
#' @inheritDotParams read_coords_tsv
#' @inheritParams tsv_paths_from_dir
#' @return A data.frame with XY
#' @export
load_tsv_masks <- function(dir_path, 
                           tsv_pattern = "^out_all_masks.tsv$",
                           position_pattern = ".*Position(\\d+)$",
                           ...){
  
  tsv_files_df <- tsv_paths_from_dir(dir_path, tsv_pattern, position_pattern)
  
  tsv_files_df.split <- split(tsv_files_df$path, tsv_files_df$pos)
  
  tsv_masks_df.split <- lapply(tsv_files_df.split, read_coords_tsv, ...)
  
  return(
    dplyr::bind_rows(tsv_masks_df.split, .id = "pos") %>% mutate(pos=as.integer(pos))
  )
}

#' Generate data.frame with paths to coordinate TSV files
#'
#' Convenience function to generate paths for \code{read_coords_tsv} and input for \code{hues_from_tsv_files2}.
#'
#' @param dir_path Path to the directory containing CellID's outputs (tipically the images' directory). CellID must have been run with TSV output enabled.
#' @param tsv_pattern A regex pattern matching the names of the TSV files.
#' @param position_pattern A regex pattern with one group for the integer position number, extracted from the directory name holding the TSV file.
#' @export
tsv_paths_from_dir <- function(dir_path,
                               tsv_pattern = "^out_all_masks.tsv$",
                               position_pattern = ".*Position(\\d+)$"){
  
  tsv_paths <- dir(unique(dir_path),
                   pattern = tsv_pattern,
                   recursive = T,
                   full.names = T)
  
  if (length(tsv_paths) == 0) {
    warning(paste("tsv_paths_from_dir: No TSV files extracted from:", dir_path))
    return(NULL)
  }
  
  tsv_dirnames <- basename(dirname(tsv_paths))
  
  tsv_files_df <-
    data.frame(
      pos = as.integer(sub(position_pattern, "\\1", tsv_dirnames)),
      path = tsv_paths
    )
  
  return(tsv_files_df)
}

#' Read TSV file with CellID mask coordinates
#' 
#' @param masks_tsv_path A path to the TSV file holding XY coordinates, from CellID's output with "-t" option.
#' @param shape_pixtype Default "b" for Hu moments based on boundary points. A character vector containing any of c("b", "i").
#' @param shape_flagtype Default 0 for Hu moments based on flag value 0. Can be any of the integer flag values present in the \code{out_bf_fl_mapping} CellID files.
#'
#' @export
read_coords_tsv <- function(masks_tsv_path, shape_pixtype = "b", shape_flagtype = 0){
  masks_coords <-readr::read_tsv(masks_tsv_path, progress = TRUE,
                                 col_types = cols(cellID = readr::col_integer(),
                                                  t.frame = readr::col_integer(),
                                                  flag = readr::col_integer(),
                                                  x = readr::col_integer(),
                                                  y = readr::col_integer(),
                                                  pixtype = readr::col_factor(levels = c("b", "i"))
                                 )
  ) %>% 
    {if(!is.null(shape_pixtype)) filter(., pixtype %in% shape_pixtype) else .} %>% 
    {if(!is.null(shape_flagtype)) filter(., flag %in% shape_flagtype) else .} %>% 
    mutate(id = factor(paste(cellID, t.frame, flag, pixtype, sep = "_")))
  
  return(masks_coords)
}

# TESTS section ####

#' Check X/Y positions of cells and cells's masks
#'
#' Make plots to examine correspondence between mean mask xpos/ypos and CellID xpos/ypos for each cell.
#'
#' @param cell_data A cell.data list object, with a "masks" dataframe element.
#'
check_tiff_mask <- function(cell_data){
  
  options(dplyr.summarise.inform = FALSE)
  
  masks_df_summary <- cell_data$masks %>% 
    group_by(cellID, t.frame, pos) %>% 
    summarise(#cellID = first(cellID),
      #t.frame = first(t.frame),
      #pos = first(pos),
      n_values = n(),
      # xspan = paste(min(x), max(x), sep = "-"),
      # yspan = paste(min(y), max(y), sep = "-"),
      xspan = min(x) - max(x) + 1,
      yspan = min(y) - max(y) + 1,
      mean_xpos = round(mean(x)),
      mean_ypos = round(mean(y))
    )
  
  cdata_compare <- left_join(mutate(cell_data$data, cellID = as.factor(cellID)),
                             select(masks_df_summary, 
                                    cellID, t.frame, pos, mean_xpos, mean_ypos), 
                             by = c("cellID", "t.frame", "pos"))
  p1 <- cdata_compare %>% 
    ggplot() +
    geom_point(aes(x = xpos, y = mean_xpos, color = factor(xpos - mean_xpos))) +
    facet_grid(pos~t.frame) +
    xlab("cellid xpos") +
    ylab("rounded mean mask xpos")
  
  p2 <- cdata_compare %>% 
    ggplot() +
    geom_point(aes(x = ypos, y = mean_ypos, color = factor(ypos - mean_ypos))) +
    facet_grid(pos~t.frame) +
    xlab("cellid ypos") +
    ylab("rounded mean mask ypos")
  
  print(p1)
  print(p2)
  
  options(dplyr.summarise.inform = TRUE)
  
  return(invisible(NULL))
}
