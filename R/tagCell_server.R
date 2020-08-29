#' Shiny app server function for tagCell
#' Define server logic required to draw a histogram
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import shiny formattable dplyr tidyr hexbin magick
#' @importFrom graphics polygon
tagCellServer <- function(input, output, session) {
  
  ### OUTPUT OBSERVERS ###
  # Reactive image 1
  output$pics <- shiny::renderImage({
    print("Rendering image 1")
    
    # Me quedo con las posiciones que me interesan
    positions <- paths$pos %>% unique()
    d <- subset(cdata, pos %in% positions & filter == TRUE)
    p <- subset(paths, pos %in% positions)

    # Output an image if filtering returns a non-empty selection
    if(nrow(d) > 0) {
      print("-- Selection not empty: magick!")
      tmpimage <- magickCell(d, p, ch=input$ch, sortVar = input$x, seed = values$seed, ...)
    } else {
      # Output white if selection is empty
      print("-- Selection is empty")
      tmpimage <- magick::image_blank(100,10,color = "white") %>% magick::image_annotate(text = "Empty set")
    }
    
    # Esto es para que shiny pueda mostrar output de imagen
    # Ver si se puede pasar a TIFF
    tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
    list(src = tmpfile, contentType = "image/jpeg")
  })
}

