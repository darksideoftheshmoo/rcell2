#' Hu moments
#' 
#' This function is CellID agnostic.
#'
#' @param xy Matriz de dos columnas (dim1 y dim2), con las coordenadas XY del boundary binario. No es importante el orden de los puntos.
#' @return Los 7 Hu moments + 1 extra (?) en un named vector
# @example
# xy <- EBImage::imageData(readImage("inst/hu_moments_img2.png"))
# xy <- which(xy>0, arr.ind = TRUE)
# plot(xy)
# xy.hu <- hu.moments(xy)
# xy.hu <- -sign(xy.hu)*log10(abs(xy.hu))
#' @export
hu.moments <- function(xy){
  ## Define function for calculating eta
  eta <- function(xy, i, j){
    ## Calculate central moment mu(i,j)
    mu <- sum(((x-cent[1])**i)*((y-cent[2])**j)*px.intensity)
    
    ## Calculate and return normalized central moment eta(i,j)
    mu/(mu00**((i+j)/2+1))
  }
  
  px.intensity <- 255
  
  x <- xy[,1]
  y <- xy[,2]
  
  ## Calculate basic object moments M(0,0) M(0,1) M(1,0) M(1,1)
  M <- matrix(sapply(1:4, function(n,i,j,x,y) sum((x**i[n])*(y**j[n])*px.intensity),
                     x=x, y=y, i=c(0,1,0,1), j=c(0,0,1,1)),2)
  
  ## Calculate object centroid (x,y)
  cent <- diag(apply(M/M[1,1],2,rev))
  
  ## Calculate object mu(0,0)
  mu00 <- sum(((x-cent[1])**0)*((y-cent[2])**0)*px.intensity)
  
  I1 <- eta(xy,2,0) + eta(xy,0,2)
  I2 <- (eta(xy,2,0) - eta(xy,0,2))**2 + 4*eta(xy,1,1)**2
  I3 <- (eta(xy,3,0) - 3*eta(xy,1,2))**2 + (3*eta(xy,2,1) - eta(xy,0,3))**2 #not very useful?
  I4 <- (eta(xy,3,0) + eta(xy,1,2))**2 + (eta(xy,2,1) + eta(xy,0,3))**2
  I5 <- (eta(xy,3,0) - 3*eta(xy,1,2))*(eta(xy,3,0) + eta(xy,1,2))*
        ((eta(xy,3,0) + eta(xy,1,2))**2 - 3*(eta(xy,2,1) + eta(xy,0,3))**2) +
        (3*eta(xy,2,1) - eta(xy,0,3)) * (eta(xy,2,1) + eta(xy,0,3)) *
        (3*(eta(xy,3,0) + eta(xy,1,2))**2 - (eta(xy,2,1) + eta(xy,0,3))**2)
  I6 <- (eta(xy,2,0) - eta(xy,0,2)) * ((eta(xy,3,0) + eta(xy,1,2))**2 -
                                       (eta(xy,2,1) + eta(xy,0,3))**2) +
        4*eta(xy,1,1) * (eta(xy,3,0) + eta(xy,1,2)) * (eta(xy,2,1) + eta(xy,0,3))
  I7 <- (3*eta(xy,2,1) - eta(xy,0,3)) * (eta(xy,3,0) + eta(xy,1,2)) *
        ((eta(xy,3,0) + eta(xy,1,2))**2 - 3*(eta(xy,2,1) + eta(xy,0,3))**2) -
        (eta(xy,3,0) - 3*eta(xy,1,2)) * (eta(xy,2,1) + eta(xy,0,3)) *
        (3*(eta(xy,3,0) + eta(xy,1,2))**2 - (eta(xy,2,1) + eta(xy,0,3))**2)
  I8 <- eta(xy,1,1) * ((eta(xy,3,0) + eta(xy,1,2))**2 - (eta(xy,0,3) + eta(xy,2,1))**2) -
        (eta(xy,2,0) - eta(xy,0,2)) * (eta(xy,3,0) + eta(xy,1,2)) * (eta(xy,0,3) + eta(xy,2,1))
  
  #hu <- c(I1, I2, I3, I4, I5, I6, I7, I8)
  #return(-sign(hu)*log10(abs(hu)))
  return(c(hum_1 = I1, 
           hum_2 = I2, 
           hum_3 = I3, 
           hum_4 = I4, 
           hum_5 = I5, 
           hum_6 = I6, 
           hum_7 = I7, 
           hum_8 = I8))
}

#' Generate a dataframe with boundary coordinates from mask tiff file
pic_df_from_tiff <- function(tiff_path, image_bits, cell_id_offset = -1){
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

#' Generate Hu moments from XY coordinates dataframe
#' 
#' This function requires a CellID column to split the coordinates by.
#' 
hues_from_xy <-  function(pic_df, split_col = "cellID"){
  # Split the dataframe by cellID
  pic_df_split <- split(pic_df, pic_df[[split_col]])
  
  # Compute Hu moments for each cellID's boundary mask XY coordinates
  hues <- lapply(pic_df_split, FUN = function(cell_coords_df){
    # Convert dataframe to matrix
    xy <- as.matrix(cell_coords_df[,c("x", "y")])
    # Rename XY columns appropriately
    colnames(xy) <- c("dim1", "dim2")
    
    # Return a named vector with the cell ids and the named hu moments
    return_vector <- c(unique(cell_coords_df[[split_col]]), rcell2::hu.moments(xy))
    names(return_vector)[1] <- split_col
    
    return_vector
  })
  
  # Bind rows from all cellIDs and return
  return(bind_rows(hues))
}

#' Generate Hu moments from XY coordinates dataframe
#' 
#' This function requires a CellID column to split the coordinates by. Optional cell-wise parallelization.
#' 
#' @param coords_df A data frame with XY coordinates, containing columns "x", "y", and the column named by the \code{split_col} parameter.
#' @param split_col A string matching the name of the colum used to split the dataframe into individual shapes.
#' 
#' @return A data frame with Hu moments for each split of the input, identified by \code{split_col}.
#' @importFrom data.table rbindlist
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @import parallel doParallel 
#' @export
#' 
hues_from_xy2 <-  function(coords_df, split_col = "cellID"){
  # Split the dataframe by cellID
  pic_df_split <- base::split(coords_df, coords_df[[split_col]], drop = T)
  split_col <- split_col
  
  # Compute Hu moments for each cellID's boundary mask XY coordinates
  n_cores <- max(1, parallel::detectCores() - 1)
  cat(paste("Message from hues_from_xy2: Parallel Hu with", n_cores, "threads.\n"))
  cl <- parallel::makeCluster(n_cores)
  # Export objects to cluster
  parallel::clusterExport(cl, "split_col", envir = environment())
  parallel::clusterExport(cl, "hu.moments")
  doParallel::registerDoParallel(cl)
  
  # Run parallelized computation
  hues <- foreach(cell_coords_df=pic_df_split) %dopar% {
    # Get "id" for the current shape
    shape_id <- unique(cell_coords_df[[split_col]])
                                
    # Convert dataframe to matrix
    xy <- as.matrix(cell_coords_df[,c("x", "y")])
    # Rename XY columns appropriately
    colnames(xy) <- c("dim1", "dim2")
    
    # Compute hu moments
    xy_humoments <- as.list(hu.moments(xy))
    xy_humoments[split_col] <- as.character.factor(shape_id)
    
    # Return as data frame
    as.data.frame(xy_humoments)
  }
  # Bind results
  hues_bound <- data.table::rbindlist(hues)
  parallel::stopCluster(cl)
  return(hues_bound)  
}

#' Read TSV file with mask coords
#' 
#' @inheritParams hues_from_tsv2
#' 
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

#' Generate Hu moments from XY coordinates TSV file
#' 
#' Cells in the dataframe are split by cellID, t.frame, flag, and pixtype by default.
#' 
#' Optional filtering by shape type (boundary or interior) and CellID flag number.
#' 
#' Optional cell-wise parallelization.
#' 
#' @param masks_tsv_path A path to the TSV file holding XY coordinates, from CellID's output with "-t" option.
#' @param .parallel Enable cell-wise parallelization using parallel::parLapply
#' @param shape_pixtype Default "b" for Hu moments based on boundary points. A character vector containing any of c("b", "i").
#' @param shape_flagtype Default 0 for Hu moments based on flag value 0. Can be any of the integer flag values present in the \code{out_bf_fl_mapping} CellID files.
#' @param cdata_subset Subset of cdata with unique rows, and only pos, t.frame and cellID columns.
#' @param position Microscope position number (int) that corresponds to current TSV file, used for filtering if \code{cdata_subset} is not NULL.
#' @import readr
#' @export
#' 
hues_from_tsv2 <- function(masks_tsv_path, .parallel = F,
                           shape_pixtype = "b", shape_flagtype = 0, 
                           cdata_subset = NULL, position = NULL){
  # Add "id" column
  cat("\nMessage from hues_from_tsv2: reading TSV id data...\n")
  masks_coords <- read_coords_tsv(masks_tsv_path, shape_pixtype, shape_flagtype)
  
  if(!is.null(cdata_subset)) 
    masks_coords <- dplyr::semi_join(mutate(masks_coords, pos = position), 
                                     cdata_subset, by = c("cellID", "pos", "t.frame"))
    
  # Compute Hu moments
  if(!.parallel) {
    hues_df <- hues_from_xy(masks_coords, split_col="id")
  } else {
    hues_df <- hues_from_xy2(coords_df = masks_coords, split_col="id")
  }
  cat(paste("Message from hues_from_tsv2: joining id data...\n"))
  hues_by_cell <- masks_coords %>% 
    select(cellID, t.frame, flag, pixtype, id) %>% unique() %>%
    left_join(hues_df, by = "id") %>%
    select(-id)
  
  cat(paste("Message from hues_from_tsv2: done!\n"))
  
  return(hues_by_cell)  # to hues_from_tsv_files2()
}

#' Generate data.frame with paths to coordinate TSV files
#' 
#' @param arguments As produced by \code{rcell2::cellArgs2}.
#' @param position_pattern A regex pattern with one group for the integer position number, extracted from the directory name holding the TSV file.
#' @param tsv_pattern A regex pattern matching the name of the TSV file.
#' @export
tsv_paths_from_args <- function(arguments,
                                position_pattern = ".*Position(\\d+)$",
                                tsv_pattern = "^out_all_masks.tsv$"){
  
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

#' Generate Hu moments data frame from one or more TSV files
#' 
#' The required TSV are generated by the mask_mod branch of the CellID program found at: https://github.com/darksideoftheshmoo/cellID-linux/tree/mask_mod
#' 
#' To play with the XY coordinates of the TIFF masks you may want to set \code{return_points = T}.
#' 
#' @param tsv_files_df A data.frame with paths to the TSV files with the cells' XY coordinates, geberated by CellID mask_mod's "-t" option. Use \code{tsv_paths_from_args()} for convenience.
#' @param return_points if TRUE it will add a "masks" dataframe to the cell_data object, containing the mask coordinates.
#' @param parralellize Enable cell-wise parallelization using parallel::parLapply
#' @param shape_pixtype Default "b" for Hu moments based on boundary points. A character vector containing any of c("b", "i").
#' @param shape_flagtype Default 0 for Hu moments based on flag value 0. Can be any of the integer flag values present in the \code{out_bf_fl_mapping} CellID files.
#' @param cdata_subset Subset of cdata with unique rows, and only pos, t.frame and cellID columns.
#' @import dplyr
#' @export
#' 
hues_from_tsv_files2 <- function(tsv_files_df, return_points = F, parralellize = T, 
                                 shape_pixtype = "b", shape_flagtype = 0, cdata_subset = NULL){
  
  if(!is.data.frame(tsv_files_df)) stop("Input error: tsv_files_df is not a data.frame.")
  if(!all(c("pos", "path") %in% names(tsv_files_df))) stop("Input error: names 'pos' and 'path' not found.")
  
  hues_df_list <- list()
  for(position in tsv_files_df[,"pos", drop = TRUE]) {
    cat("\nComputing Hu moments for position:", position)
    
    masks_tsv_path <- tsv_files_df %>% 
      filter(pos == position)%>% 
      .[,"path", drop = TRUE] %>% 
      {
        if(length(.) != 1) stop("Error in append_hues2: more than one TSV file per position") 
        else .
      }
  
    hues_df <- 
      hues_from_tsv2(masks_tsv_path = masks_tsv_path,
                     .parallel = parralellize, 
                     cdata_subset = cdata_subset,
                     shape_pixtype = shape_pixtype, 
                     shape_flagtype = shape_flagtype,
                     position = position) %>% 
      mutate(pos = as.integer(position))
    
    hues_df_list[[position]] <- hues_df
  }
  
  return(bind_rows(hues_df_list))
}


#' Add Hu moments data.frame to a cell_data object
#' 
#' The required TSV are generated by the mask_mod branch of the CellID program found at: https://github.com/darksideoftheshmoo/cellID-linux/tree/mask_mod
#' 
#' To play with the XY coordinates of the TIFF masks you may want to set \code{return_points = T}.
#' 
#' @param tsv_files_df A data.frame with paths to the TSV files with the cells' XY coordinates, geberated by CellID mask_mod's "-t" option.
#' @param cell_data NULL list by default. You may pass a cell.data list object from Rcell2 load_cell_data, for convenience, or just the cdata dataframe as \code{list(data = cdata)}.
#' @param return_points if TRUE it will add a "masks" dataframe to the cell_data object, containing the mask coordinates.
#' @param shape_pixtype Default "b" for Hu moments based on boundary points. A character vector containing any of c("b", "i").
#' @param shape_flagtype Default 0 for Hu moments based on flag value 0. Can be any of the integer flag values present in the \code{out_bf_fl_mapping} CellID files.
#' @param overwrite Default FALSE, overwrite Hu moments element in the input list if already present.
#' @export
#' 
#' @return The Hu moments dataframe usis assgned to the input list as an element named "Hu_moments".
#' 
append_hues2 <- function(tsv_files_df, 
                         cell_data = NULL, 
                         return_points = F,
                         parralellize = T,
                         shape_pixtype = "b",
                         shape_flagtype = 0,
                         overwrite = F
                         ){
  
  if(return_points) 
    stop("return_points option not implemented yet :/")
  if("Hu_moments" %in% names(cell_data) && !overwrite) 
    stop("cell_data already has a Hu_moments element, use 'overwrite=TRUE' to force it.")
  
  if(!is.null(cell_data)) 
    cdata_subset <- unique(cell_data[["data"]][, c("cellID", "pos", "t.frame")]) 
  else 
    cdata_subset <- NULL
  
  hues_df <- hues_from_tsv_files2(tsv_files_df, 
                                  return_points = F, 
                                  parralellize = parralellize,
                                  shape_pixtype = shape_pixtype, 
                                  shape_flagtype = shape_flagtype, 
                                  cdata_subset = cdata_subset)
  
  # Append Hu moment data to cell_data object
  if(is.null(cell_data)) list(Hu_moments = hues_df) 
  else cell_data[["Hu_moments"]] <- hues_df

  # if(return_points){
  #   masks_df <- bind_rows(
  #     lapply(pic.and.hues.dfs, FUN = function(l) l[["pic_df"]])
  #   )
  # 
  #   cell_data[["masks"]] <- masks_df
  # }

  return(cell_data)
}

#' Append Hu moments columns to data on a cell_data object
#' 
#' It is imperative that the TIFF images contain only cell boundaries, with pixel intensity values such that "CellID = (2^image_bits -1) - pixel_value"
#' 
#' This is automatically generated by the mask_mod branch of the CellID program found at: https://github.com/darksideoftheshmoo/cellID-linux/tree/mask_mod
#' 
#' To play with the XY coordinates of the TIFF masks you may want to set \code{return_points = T}.
#' Or if you want to skip calculation of the Hu moments, checkout \code{rcell2:::pic_df_from_tiff}.
#' 
#' @param cell_data a cell_data object as loaded by rcell2::load_cell_data
#' @param image_bits an integer indicating the bit-depth on the TIFF images, such that maximum intensity is \code{image_bits^2 -1} and minimum is zero.
#' @param cell_id_offset the offset respect to maximum pixel intensity, such that \code{cellID = maximum_intensity - boundary_intensity + cell_id_offset}.
#' @param return_points if TRUE it will add a "masks" dataframe to the cell_data object, containing the mask coordinates.
#' @export
#' 
append_hues <- function(cell_data, image_bits, cell_id_offset = -1, return_points = F){
  paths <- cell_data$images %>% 
    mutate(file = paste0(path, "/", image)) %>% 
    filter(channel == "BF.out")
  
  # For testing
  # pic_metadata <- c()
  # pic_metadata["file"] <- paths[1,]$file
  # pic_metadata["t.frame"] <- paths[1,]$t.frame
  # pic_metadata["pos"] <- paths[1,]$pos
  # image_bits = 16
  # return_points = T
  
  pic.and.hues.dfs <- apply(paths, MARGIN = 1, FUN = function(pic_metadata){
    # Extract xy coordinates list from the BF.out mask tiff and derive cellID from the intensities
    pic_df <- pic_df_from_tiff(tiff_path = pic_metadata["file"], 
                               image_bits = image_bits, 
                               cell_id_offset = cell_id_offset)
    
    # Compute Hu moments for each cellID
    hues_df <- hues_from_xy(pic_df)
    
    # Add position and time information
    hues_df$t.frame <- as.integer(pic_metadata["t.frame"])
    hues_df$pos <- as.integer(pic_metadata["pos"])
    pic_df$t.frame <- as.integer(pic_metadata["t.frame"])
    pic_df$pos <- as.integer(pic_metadata["pos"])
    
    # Return the Hu moments and points dataframes
    list(pic_df = if(return_points) pic_df else NULL, 
         hues_df = hues_df)
  })
  
  # Bind hues dataframes from all pics
  hues.df <- bind_rows(
    lapply(pic.and.hues.dfs, FUN = function(l) l[["hues_df"]])
  )
  
  # Append Hu moment data to cell_data object
  cell_data[["data"]] <- left_join(cell_data[["data"]],
                                   hues.df, 
                                   by = c("cellID", "t.frame", "pos"))
  
  if(return_points){
    masks_df <- bind_rows(
      lapply(pic.and.hues.dfs, FUN = function(l) l[["pic_df"]])
    )
    
    cell_data[["masks"]] <- masks_df
  }
  
  return(cell_data)
}

#' Make plots to examine correspondence between mean mask xpos/ypos and CellID xpos/ypos for each cell.
check_tiff_mask <- function(cell_data){
  
  options(dplyr.summarise.inform = FALSE)
  
  masks_df_summary <- cell_data$masks %>% 
    group_by(cellID, t.frame, pos) %>% 
    summarise(cellID = first(cellID),
              t.frame = first(t.frame),
              pos = first(pos),
              n_values = n(),
              # xspan = paste(min(x), max(x), sep = "-"),
              # yspan = paste(min(y), max(y), sep = "-"),
              xspan = min(x) - max(x),
              yspan = min(y) - max(y),
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

#' Make plots to examine correspondence between mean mask xpos/ypos and CellID xpos/ypos for each cell.
check_tiff_mask2 <- function(cell_data){
  
  options(dplyr.summarise.inform = FALSE)
  
  masks_df_summary <- cell_data$masks %>% 
    group_by(cellID, t.frame, pos) %>% 
    summarise(cellID = first(cellID),
              t.frame = first(t.frame),
              pos = first(pos),
              n_values = n(),
              # xspan = paste(min(x), max(x), sep = "-"),
              # yspan = paste(min(y), max(y), sep = "-"),
              xspan = min(x) - max(x),
              yspan = min(y) - max(y),
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