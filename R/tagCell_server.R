#' Shiny app server function for tagCell
#' Define server logic required to draw a histogram
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import shiny formattable dplyr tidyr hexbin magick
#' @importFrom graphics polygon
tagCellServer <- function(input, output, session) {
  
  print("Appending tags to tempfile:")
  print(tmp_csv_output)
  
  d <- cdata %>% dplyr::arrange(ucid)
  p <- paths
  
  reactive_values <- shiny::reactiveValues(ith_cell = 1,
                                           i_line = 1,
                                           other_reactive_values = c())
  
  ### BUTTON OBSERVERS ###
  shiny::observeEvent(
    eventExpr = input$next_cell,
    handlerExpr = {
      reactive_values$ith_cell = min(length(unique(cdata$ucid)), reactive_values$ith_cell + 1)
    })
  
  shiny::observeEvent(
    eventExpr = input$prev_cell,
    handlerExpr = {
      reactive_values$ith_cell = max(1, reactive_values$ith_cell - 1)
    })
  
  # BUTTON OBSERVER 2: EXIT  ----------------
  observeEvent(
    # Acción al apretar el botón de cerrar la app
    eventExpr = input$quit,
    handlerExpr = {
      writeLines("\nQuit event fired!")
      stopApp(tmp_csv_output)
    }
  )
  
  ### OUTPUT OBSERVERS and RENDERERS ###
  # Reactive text 1
  output$cell_ith <- shiny::renderText({
    ith_cell_ucid <- d$ucid[reactive_values$ith_cell]
    
    tag_line <- paste(
      shiny::isolate(reactive_values$i_line),
      ith_cell_ucid,
      sep = ", "
    )
    
    write(tag_line, file=tmp_csv_output,append=TRUE)
    
    shiny::isolate({reactive_values$i_line <- reactive_values$i_line +1})
    
    reactive_values$ith_cell
  })
  
  # Reactive image 1
  output$pics <- shiny::renderImage({
    print("Rendering image 1")
    
    if(nrow(d) > 0) {
      print("-- Selection not empty: magick!")
      cdata.selected <- d[d$ucid == d$ucid[reactive_values$ith_cell],]
      magick.cell <-  magickCell(cdata = cdata.selected, 
                                 p,
                                 ch=input$image_channel, 
                                 n = n_max, 
                                 .equalize = F,
                                 .normalize = T)
      tmpimage <- magick.cell$img
      print(magick.cell$ucids)
    } else {
      # Output white if selection is empty
      print("-- Selection is empty")
      tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
    }
    
    tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
    list(src = tmpfile, contentType = "image/jpeg")
  }, deleteFile=TRUE)
  
  # Reactive plot 1
  output$plot <- shiny::renderPlot({
    print("Rendering image 1")
    
    if(!is.null(tag_ggplot)){
      tag_ggplot
    }
  })
}


