#' Extract additional info from a `TIFFTAG_DESCRIPTION`.
#'
#' @return A list with the `description` attribute modified.
#'
#' @noRd
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
      descr <- matrix(unlist(strsplit(y[1,2],"(&amp;#13;&amp;#10;)|(: )")),ncol=2,byrow=TRUE)
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
#' TIFF files contain metadata about images in their _TIFF tags_. This function
#' is for reading this information without reading the actual image. This is an
#' extention of the \code{read_tags} function from the \code{ijtiff} package
#' for extracting additional information nested within the _description_
#' attribute of some TIFF images.
#'
#' @param path A string. The path to the tiff file to read.
#' @param frames Which frames do you want to read tags from. Default first frame
#'   only. To read from the 2nd and 7th frames, use `frames = c(2, 7)`, to read
#'   from all frames, use `frames = "all"`.
#'
#' @return A list of lists.
#'
#' @author Simon Urbanek, Kent Johnson, Rory Nolan, Andreas Constantinou.
#'
#' @examples
#' rcell2::read_tiff_tags(path = "data/image_samples/BF_Position001.tif", frames = 1)
#' 
#' @export
read_tiff_tags <- function(path,frames=1){
  lapply(rcell2::read_tags(path,frames),read_description)
}
