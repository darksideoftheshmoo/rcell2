read_bf_masks <- function(path, cell_id_offset = -1){
  # Read image
  img <- ijtiff::read_tif(path = path,msg = FALSE)
  
  # Get image bitsize and get max image cellID intensity value
  image_bits <- attributes(img)$bits_per_sample
  stopifnot("Expected 16-bit BF.out image." = image_bits==16)
  max_intensity <- 2**image_bits-1+cell_id_offset
  
  # Get coordinates and values for all pixels with intensities > 0
  mask_coord <- which(img>0,arr.ind=TRUE)
  mask_vals <- img[mask_coord]
  
  # Get lowest intensity value
  min_intensity <- min(mask_vals)
  
  # Calculate center intensity
  center_intensity <- (max_intensity)-(max_intensity-min_intensity)%/%2
  
  # Check if center_intensity exists in image
  if(is.na(match(center_intensity,mask_vals))){
    # If center_intensity does not exist, then there are interior pixels offset
    # from boundary pixels
    
    # Get boundary indices, coordinates, and values
    boundary_idx <- which(mask_vals>center_intensity)
    boundary_coord <- mask_coord[boundary_idx,1:2]
    boundary_vals <- mask_vals[boundary_idx]
    
    # Get interior indices, coordinates, and values
    interior_idx <- which(mask_vals<center_intensity)
    interior_coord <- mask_coord[interior_idx,1:2]
    interior_vals <- mask_vals[interior_idx]
    interior_shift <- max_intensity-max(interior_vals)
    
    # Make mask_data matrix and insert boundary and interior data
    mask_data <- matrix(0,nrow=nrow(mask_coord),ncol=5)
    mask_data[,1] <- c(-boundary_vals,-(interior_vals+interior_shift))+max_intensity
    mask_data[,2:3] <- rbind(boundary_coord,interior_coord)
    mask_data[,4] <- c(boundary_vals,interior_vals)
    
    # Set pixel type (boundary = 1, interior = 2)
    mask_data[,5] <- c(rep(1,nrow(boundary_coord)),rep(2,nrow(interior_coord)))
    
  } else {
    # If center_intensity does exists, then there are no interior pixels offset
    # from boundary pixels
    
    # Make mask_data matrix and insert boundary and interior data
    # Pixel type defaults to 0
    mask_data <- matrix(0,nrow=nrow(mask_coord),ncol=5)
    mask_data[,1] <- max_intensity-mask_vals
    mask_data[,2:3] <- mask_coord[,1:2]
    mask_data[,4] <- mask_vals
  }
  # Set column names, convert to data.frame and return
  mask_data <- mask_data[order(mask_data[,1]),]
  colnames(mask_data) <- c("cellID","x","y","pix_value","pix_type")
  as.data.frame(mask_data)
}

img_path <- "~/Mikroskopi/201025_example_dataset2/no_interior_offset/BF_Position001_time001.tif.out.tif"
test1 <- read_bf_masks(path = img_path)

img_path <- "~/Mikroskopi/201025_example_dataset2/with_interior_offset/BF_Position001_time001.tif.out.tif"
test2 <- read_bf_masks(path = img_path)
