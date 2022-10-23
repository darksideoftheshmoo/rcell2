#' Extract additional info from a `TIFFTAG_DESCRIPTION`.
#'
#' @return A list with the `description` attribute modified.
#'
#' @keywords internal
read_description <- function(x){
  # If "description" attribute contains HTML-formatted metadata
  if(length(grep("<MetaData>",x$description))){
    # Split out properties (<prop ) and custom properties (<custom-prop )
    y <- unlist(strsplit(x$description,"(<prop )|(<custom-prop )|(</prop>)|(</custom-prop>)"))
    
    # Filter for rows with "id=" that contain properties
    y <- y[grep("(id=)",y)]
    
    # Split each property at "id=", " type=", and " value="
    y <- unlist(strsplit(y,"(id=)|( type=)|( value=)"))
    
    # Clean rows of special characters and convert to matrix
    y <- matrix(gsub("(>\r\n\t\t)|(>\r\n\t)|(\")","",y),ncol=4,byrow=TRUE)[,-c(1,3)]
    
    # Check if first row ("Description") contains HTML character codes
    if(length(grep("(&amp;#13;&amp;#10;)",y[1,2]))){
      # Extract and split value in first row
      descr <- matrix(data=unlist(strsplit(y[1,2],"(&amp;#13;&amp;#10;)|(: )")),
                      ncol=2,
                      byrow=TRUE)
      y <- rbind(descr,y[-1,])
    }
    
    # Convert to data.frame, set column names, add to tags list
    y <- as.data.frame(y)
    names(y) <- c("variable","value")
    x$description <- y
  }
  # Else if "description" attribute contains newline (\n) characters
  else if(length(grep("\n",x$description))){
    # Split out rows and columns
    y <- matrix(unlist(strsplit(x$description,"(\n)|(=)")),ncol=2,byrow=TRUE)
    y <- as.data.frame(y)
    names(y) <- c("variable","value")
    x$description <- y
  }
  x
}

#' Read TIFF tag information.
#' 
#' TIFF files contain metadata about images in their \emph{TIFF tags}. This function
#' is for reading this information without reading the actual image. It extends
#' on the \code{\link[ijtiff]{read_tags}} function from the \link{ijtiff} package, and extracts
#' additional information nested within the \code{description} attribute of some TIFF
#' images.
#'
#' @param path A string. The path to the tiff file to read.
#' @param frames Which frames do you want to read tags from. Default first frame
#'   only. To read from the 2nd and 7th frames, use `frames = c(2, 7)`, to read
#'   from all frames, use `frames = "all"`.
#'
# @importFrom ijtiff read_tags
#'
#' @return A list of lists.
#'
#' @examples
#' rcell2::read_tiff_tags(path = "data/image_samples/BF_Position001.tif", frames = 1)
#' 
#' @export
read_tiff_tags <- function(path,frames=1){
  # Check dependency  
  if(!"ijtiff" %in% rownames(installed.packages())){
    stop("read_tiff_tags requires the 'ijtiff' package, which is not installed by default. Please install it to use this function :)")
  }
  lapply(ijtiff::read_tags(path,frames), read_description)
  # stop("read_tiff_tags: the dependency ijtiff::read_tags was removed due to compilation errors.")
}

#' Get stage XY coordinates from tiff file's description tag
#' 
#' The description is in the tiff file's metadata, tipically written by Metamorph.
#' 
#' @inheritParams get_tiff_description
#' @return A named vector, with "x" and "y" coordinates.
#' @export
#' @examples 
#' 
#' # Get the "images" data.frame from CellID
#' cell.data <- cell.load.alt("path/to/output/")
#' images <- cell.data$images
#' 
#' # Generate a list with file paths split and named by "pos"
#' tags.list <- images %>% filter(t.frame==0, channel=="BF") %>% 
#'   {split(.$file, .$pos)}
#'   
#' # Otherwise, a simple list of file paths ordered by "pos" suffice.
#' # tags.list <- list(file.paths)  # NOT RUN
#' 
#' # Get stage positions
#' tags.stage <- tags.list %>% 
#'   lapply(get_stage_xy_from_file) %>% bind_rows(.id="pos") %>% 
#'   mutate(pos = as.integer(pos))
#' 
#' tags.stage
#' 
get_stage_xy_from_file <- function(tiff.path, frames=1){
  
  tags.xml <- get_tiff_description(tiff.path, frames)
  
  tag.attrs <- tags.xml$MetaData$PlaneInfo %>% lapply(attributes) %>% 
    {setNames(., sapply(., "[[", "id"))}
  
  c(
    x = tag.attrs$`stage-position-x`$value %>% as.numeric(),
    y = tag.attrs$`stage-position-y`$value %>% as.numeric()
  )
}

#' Get tiff acquisition times from tiff file's description tag
#' 
#' The description is in the tiff file's metadata, tipically written by Metamorph.
#' 
#' @inheritParams get_tiff_description 
#' @inheritParams base::as.POSIXct
#' @param ... Passed to \code{as.POSIXct}.
#' @export
#' @return The tiff time (as.POSIXct).
#' 
get_tiff_time_from_file <- function(tiff.paths, frames=1,
                                    format = "%Y%m%d %H:%M:%OS",
                                    ...){
  tiff.times <- 
    sapply(tiff.paths, FUN = function(tiff.path){
      tags.xml <- get_tiff_description(tiff.path, frames)
      
      tag.attrs <- tags.xml$MetaData$PlaneInfo %>% lapply(attributes) %>% 
        {setNames(., sapply(., "[[", "id"))}
      
      time.str <- tag.attrs$`acquisition-time-local`$value
      
      return(time.str)
    })
  
  tiff.times.df <- data.frame(time.str = tiff.times) %>% 
    mutate(time.posixct = as.POSIXct(x = time.str,
                                     format=format,
                                     ...))
  
  return(tiff.times.df)
}

#' Get Metamorph image XML description
#' 
#' @param tiff.path A path to the tiff file, passed to \code{\link[ijtiff]{read_tags}}.
#' @param frames The frame number to process, passed to \code{\link[ijtiff]{read_tags}}.
#' 
#' 
get_tiff_description <- function(tiff.path, frames=1){
  
  # Check dependency  
  if(!requireNamespace("ijtiff")){
    stop("get_tiff_description requires the 'ijtiff' package.")
  }
  
  description <- tiff.path %>% 
    ijtiff::read_tags(frames = frames) %>% 
    .[["frame1"]] %>% .[["description"]] %>% 
    xml2::read_xml() %>% xml2::as_list()
  
  return(description)
}
