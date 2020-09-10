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
  write("", file=tmp_csv_output,append=TRUE)
  
  d <- cdata %>% dplyr::arrange(ucid, t.frame) %>% 
    mutate(ucid_t.frame = paste(ucid, t.frame, sep = "_"))
  p <- paths
  
  reactive_values <- shiny::reactiveValues(ith_cell = 1,
                                           i_line = 1,
                                           other_reactive_values = c(),
                                           selected_cell_tags = list())
  
  ### UI OBSERVERS   ----------------
  output$moreControls <- renderUI({
    # tagList(
      # sliderInput("n", "N", 1, 1000, 500),
      # textInput("label", "Label")
    # )
    # shiny::selectInput('tags','Tag', cell_tags, multiple = TRUE, selected = NULL, selectize = T)
    # ith_cell <- reactive_values$ith_cell
    # selected_cell_tags <- reactive_values$selected_cell_tags
    # 
    # print("- Generating seelctInput fields...")
    # print(ith_cell)
    # print(selected_cell_tags[ith_cell])
    
    lapply(1:length(names(cell_tags)), function(tag_group){
      # shiny::selectInput(inputId = names(cell_tags)[tag_group], 
      #                    label = names(cell_tags)[tag_group],
      #                    choices = cell_tags[tag_group],
      #                    multiple = TRUE,
      #                    # selected = unlist(selected_cell_tags[ith_cell])[tag_group],
      #                    selected = NULL,
      #                    selectize = T)
      shiny::checkboxGroupInput(inputId = names(cell_tags)[tag_group], 
                                label = names(cell_tags)[tag_group],
                                choices = unlist(cell_tags[tag_group])
                                # choices = list("asd", "dddd")
                                )
    })
  })
  
  ### BUTTON OBSERVERS   ----------------
  # BUTTON 1.1: NEXT  ----------------
  shiny::observeEvent(
    eventExpr = input$next_cell,
    handlerExpr = {
      print("- Next cell requested, saving current tags...")
      ith_cell <- reactive_values$ith_cell                         # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid for that cell
      
      print("-- Saving tag selection")
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
          ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      # Handle previous > total
      reactive_values$ith_cell <- ith_cell + 1                     # Update the ith_cell reactive value
      if(ith_cell + 1 > nrow(cdata)){
        showNotification("There is no next cell, staying at the first one.")
        reactive_values$ith_cell <- 1
      }
    })
  
  # BUTTON 1.2: PREVIOUS  ----------------
  shiny::observeEvent(
    eventExpr = input$prev_cell,
    handlerExpr = {
      print("- Previous cell requested...")
      ith_cell <- reactive_values$ith_cell                       # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell]) # Get ucid for that cell
      
      print("-- Saving tag selection")
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
          ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      # Handle previous < 1
      reactive_values$ith_cell <- ith_cell - 1                   # Update the ith_cell reactive value
      if(reactive_values$ith_cell < 1){
        showNotification("There is no previous cell, staying at the first one.")
        reactive_values$ith_cell <- 1
      }
    })
  
  # SIDE EFFECTS FOR PREV/NEXT BUTTON
  shiny::observe({
    print("-- Updating tag selection for next cell")
    selected_cell_tags <- reactive_values$selected_cell_tags
    ith_cell <- reactive_values$ith_cell
    ith_ucid <- as.character(d$ucid_t.frame[ith_cell])
    
    if(ith_ucid %in% names(selected_cell_tags)){
      print("--- UCID tag found")
      # selected_cell_tags[[ith_ucid]]
      selected_cell_tags <- selected_cell_tags[[ith_ucid]]
      for(tag_group in names(cell_tags)){
        if(tag_group %in% names(selected_cell_tags)){
          # shiny::updateSelectInput(session,
          shiny::updateCheckboxGroupInput(session,
                                          inputId = tag_group,
                                          choices = cell_tags[[tag_group]],
                                          selected = selected_cell_tags[[tag_group]])
        } else {
          # shiny::updateSelectInput(session,
          shiny::updateCheckboxGroupInput(session,
                                          inputId = tag_group,
                                          choices = cell_tags[[tag_group]],
                                          selected = NULL)
        }
      }
    } else {
      print("--- UCID not tagged")
      for(tag_group in names(cell_tags)){
        # shiny::updateSelectInput(session,
        # shiny::updateSelectInput(session,
        shiny::updateCheckboxGroupInput(session,
                                        inputId = tag_group,
                                        choices = cell_tags[[tag_group]],
                                        selected = NULL)
      }
    }
  })
  
  # BUTTON OBSERVER 2: EXIT  ----------------
  observeEvent(
    # Acción al apretar el botón de cerrar la app
    eventExpr = input$quit,
    handlerExpr = {
      writeLines("\nQuit event fired!")
      
      output <- reactive_values$selected_cell_tags %>% 
        bind_rows(.id = "ucid_t.frame") %>% #%>% mutate(ucid = as.numeric(ucid_t.frame))
        separate(ucid_t.frame, c("ucid", "t.frame"))

      
      # stopApp(list(tmp_csv_output, reactive_values$selected_cell_tags))
      # stopApp(tmp_csv_output)
      stopApp(output)
    }
  )
  
  ### OUTPUT OBSERVERS and RENDERERS  ----------------
  # Reactive text 1  ----------------
  output$cell_ith <- shiny::renderText({
    ith_ucid <- d$ucid[reactive_values$ith_cell]
    
    tag_line <- paste(
      shiny::isolate(reactive_values$ith_cell),
      ith_ucid,
      # paste(unlist(reactive_values$selected_cell_tags[[ith_ucid]]), collapse = "; "),
      sep = ", "
    )
    
    write(tag_line, file=tmp_csv_output,append=TRUE)
    
    paste(
    "ucid: ", d[reactive_values$ith_cell, c("ucid")],
    "\nt.frame: ", d[reactive_values$ith_cell, c("t.frame")]
    )
  })
  
  # Reactive image 1: magickCell  ----------------
  output$pics <- shiny::renderImage({
    print("- Rendering image 1")
    
    if(nrow(d) > 0) {
      print("-- Selection not empty: magick!")
      # cdata.selected <- d[d$ucid == d$ucid[reactive_values$ith_cell],]
      cdata.selected <- d[reactive_values$ith_cell,]
      magick.cell <-  magickCell(cdata = cdata.selected, 
                                 p,
                                 # ch=input$image_channel, 
                                 cell_resize=cell_resize,
                                 ch=tag_channels_select, 
                                 n = n_max, 
                                 .equalize = F,
                                 .normalize = T,
                                 boxSize = tag_box_size, 
                                 return_single_imgs = T)
      tmpimage <- magick.cell$img
      print(magick.cell$ucids)
    } else {
      # Output white if selection is empty
      print("-- Selection is empty")
      tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
    }
    
    tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
    
    list(src = tmpfile, 
         contentType = "image/jpeg")
  }, deleteFile=TRUE)
  
  # Reactive plot 1  ----------------
  output$plot <- shiny::renderPlot({
    print("- Rendering plot 1")
    
    ith_ucid <- as.character(d$ucid[reactive_values$ith_cell])
    ith_t.frame <- as.character(d$t.frame[reactive_values$ith_cell])
    
    print(ith_ucid)
    
    ucid_data <- filter(cdata, ucid == ith_ucid)
    
    if(!is.null(tag_ggplot)){
      tag_ggplot <- tag_ggplot %+% ucid_data
      tag_ggplot + geom_vline(xintercept = as.numeric(ith_t.frame), color = "red")
    }
  })
}


