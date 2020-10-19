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
#' This function requires a CellID column to split the coordinates by.
#' 
hues_from_xy2 <-  function(pic_df, split_col = "cellID"){
  # Split the dataframe by cellID
  pic_df_split <- split(pic_df, pic_df[[split_col]])
  
  # Compute Hu moments for each cellID's boundary mask XY coordinates
  cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
  hues <- parallel::parLapply(cl, 
                              pic_df_split, 
                              fun = function(cell_coords_df){
    # Convert dataframe to matrix
    xy <- as.matrix(cell_coords_df[,c("x", "y")])
    # Rename XY columns appropriately
    colnames(xy) <- c("dim1", "dim2")
    
    # Return a named vector with the cell ids and the named hu moments
    return_vector <- c(unique(cell_coords_df[[split_col]]), rcell2::hu.moments(xy))
    names(return_vector)[1] <- split_col
    
    return_vector
  })
  parallel::stopCluster(cl)
  
  # Bind rows from all cellIDs and return
  return(bind_rows(hues))
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