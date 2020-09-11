#' Shiny app server function for tagCell
#' Define server logic required to draw a histogram
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import shiny formattable dplyr tidyr hexbin magick
#' @importFrom graphics polygon
tagCellServer <- function(input, output, session) {
  
  if(is.null(tmp_output_file)){
    tmp_output_file <- tempfile(tmpdir = "./", fileext = ".txt")
  }
  print(paste("Appending tags to tempfile:", tmp_output_file))
  dir.create(dirname(normalizePath(tmp_output_file)), recursive = T)
  # write("", file=tmp_output_file)
  
  d <- cdata %>% dplyr::arrange(ucid, t.frame) %>% 
    mutate(ucid_t.frame = paste(ucid, t.frame, sep = "_")) %>% 
    mutate(
      cellID = as.integer(cellID),
      ucid = as.integer(ucid),
      t.frame = as.integer(t.frame)
    )
  p <- paths
  
  ucid.unique <- unique(d$ucid)
  
  reactive_values <- shiny::reactiveValues(ith_cell = 1, ith_ucid = numeric(),
                                           i_line = 1,
                                           other_reactive_values = c(),
                                           selected_cell_tags = list(),
                                           click_vars = list())
  
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
  
  ### INPUT OBSERVERS   ----------------
  # PLOT CLICK OBSERVER   ----------------
  observeEvent(input$plot_click, handlerExpr = {
    print("- Plot clicked")
    .variables <- list(x = isolate(input$vertex1$x),
                       y = isolate(input$vertex1$y),
                       xvar = quo_name(tag_ggplot$layers[[1]]$mapping$x),
                       yvar = quo_name(tag_ggplot$layers[[1]]$mapping$y))
    
    reactive_values$click_vars <- .variables
  })
  # BUTTON 1.1: NEXT  ----------------
  shiny::observeEvent(
    eventExpr = input$next_cell,
    handlerExpr = {
      print("- Next cell requested, saving current tags...")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_cell")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid for that cell
      
      print("-- Saving tag selection")
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
          ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      # reactive_values$ith_ucid <- as.character(d$ucid[reactive_values$ith_cell])  # Save the ucid
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags             # Save the tag list to a UCID name element in a reactive values list.
      
      # Handle previous > total
      if(ith_cell == nrow(cdata)){
        showNotification("There is no next cell, staying at the current one.", type = "warning")
      } else {
        reactive_values$ith_cell <- ith_cell + 1                     # Update the ith_cell reactive value
      }
    })
  
  # BUTTON 1.2: PREVIOUS  ----------------
  shiny::observeEvent(
    eventExpr = input$prev_cell,
    handlerExpr = {
      print("- Previous cell requested...")
      shinyjs::disable("prev_cell")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                               # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell]) # Get ucid for that cell
      
      print("-- Saving tag selection")
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
          ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      # Handle previous < 1
      if(reactive_values$ith_cell < 1){
        showNotification("There is no previous cell, staying at the current one.", type = "warning")
      } else {
        reactive_values$ith_cell <- ith_cell - 1                   # Update the ith_cell reactive value
      }
    })
  
  # BUTTON 2.1: NEXT UCID  ----------------
  shiny::observeEvent(
    eventExpr = input$next_ucid,
    handlerExpr = {
      print("- Next ucid requested, saving current tags...")
      shinyjs::disable("next_ucid")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_cell")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid_t.frame for that cell
      
      print("-- Saving tag selection")
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
        ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      # Skip to the next UCID
      ucid.oi <- as.character(d$ucid[reactive_values$ith_cell])   # Get bare ucid for the current cell
      ucid.oi.index <- match(ucid.oi, ucid.unique)
      ucid.next <- ucid.unique[ucid.oi.index + 1]                 # Get the next ucid
      ucid.next.index <- match(ucid.next, d$ucid)                 # Get the next ucid's row index
      # Handle next > total
      if(ucid.oi.index >= length(ucid.unique)){
        showNotification("--- There is no next ucid, staying at the current one.", type = "warning")
      } else {
        showNotification("--- Moving to next ucid.", duration = 1)
        reactive_values$ith_cell <- ucid.next.index               # Update the ith_cell reactive value
      }
    })
  
  # BUTTON 2.2: PREVIOUS UCID ----------------
  shiny::observeEvent(
    eventExpr = input$prev_ucid,
    handlerExpr = {
      print("- Previous ucid requested...")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("prev_cell")
      shinyjs::disable("next_cell")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                               # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell]) # Get ucid for that cell
      
      print("-- Saving tag selection")
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
        ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      # Skip to the previous UCID
      ucid.oi <- as.character(d$ucid[reactive_values$ith_cell])   # Get bare ucid for the current cell
      ucid.oi.index <- match(ucid.oi, ucid.unique)
      ucid.next <- ucid.unique[ucid.oi.index - 1]                 # Get the previous ucid
      ucid.next.index <- match(ucid.next, d$ucid)                 # Get the previous ucid's row index
      
      # Handle next > total
      if(ucid.oi.index == 1){
        showNotification("--- There is no previous ucid, staying at the current one.", type = "warning")
      } else {
        showNotification("--- Moving to previous ucid.", duration = 1)
        reactive_values$ith_cell <- ucid.next.index               # Update the ith_cell reactive value
      }
    })
  
  # SIDE EFFECTS FOR PREV/NEXT BUTTON  ----------------
  shiny::observe({
    print("-- Updating tag selection for next or previous cell")
    shinyjs::disable("prev_cell")
    shinyjs::disable("next_cell")
    shinyjs::disable("prev_ucid")
    shinyjs::disable("next_ucid")
    shinyjs::disable("moreControls")
    
    ith_cell <- reactive_values$ith_cell
    ith_ucid <- as.character(d$ucid_t.frame[ith_cell])
    selected_cell_tags <- isolate(reactive_values$selected_cell_tags)
    
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
      print("--- UCID not yet tagged")
      for(tag_group in names(cell_tags)){
        # shiny::updateSelectInput(session,
        # shiny::updateSelectInput(session,
        shiny::updateCheckboxGroupInput(session,
                                        inputId = tag_group,
                                        choices = cell_tags[[tag_group]],
                                        selected = NULL)
      }
    }
    
    shinyjs::delay(300, expr = {
      shinyjs::enable("prev_cell")
      shinyjs::enable("next_cell")
      shinyjs::enable("prev_ucid")
      shinyjs::enable("next_ucid")
      shinyjs::enable("moreControls")
      
    })
  })
  
  # BUTTON OBSERVER 3: EXIT  ----------------
  observeEvent(
    # Acción al apretar el botón de cerrar la app
    eventExpr = input$quit,
    handlerExpr = {
      writeLines("\n- Quit event fired")

      print(paste("-- Saving progress to file:", tmp_output_file))
      
      table_output <- reactive_values$selected_cell_tags %>% bind_rows(.id = "ucid_t.frame")
      if(nrow(table_output) > 0){
        table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
          mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
          left_join(select(d, ucid, pos, cellID))
        table_output %>% readr::write_csv(path = tmp_output_file)
      } else {
        table_output <- data.frame()
      }
      
      print("-- Returning progress to output:")
      # stopApp(list(tmp_csv_output, reactive_values$selected_cell_tags))
      # stopApp(tmp_csv_output)
      stopApp(table_output)
    }
  )
  
  # BUTTON OBSERVER 4: SAVE  ----------------
  observeEvent(
    # Acción al apretar el botón de cerrar la app
    eventExpr = input$save,
    handlerExpr = {
      writeLines("\n- Save event fired")
      
      table_output <- reactive_values$selected_cell_tags %>% 
        bind_rows(.id = "ucid_t.frame")  #%>% mutate(ucid = as.numeric(ucid_t.frame))
      
      if(nrow(table_output) > 0){
        table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
          mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
          left_join(select(d, ucid, pos, cellID))
        showNotification(paste("-- Saving progress to file:", tmp_output_file), duration = 4, type = "message")
      } else {
        table_output <- data.frame(message = "No annotations yet...")
        showNotification(paste("-- No annotations yet, nothing was saved."), duration = 4, type = "message")
      }
      
      readr::write_csv(table_output, path = tmp_output_file)
      
    }
  )
  
  ### OUTPUT OBSERVERS and RENDERERS  ----------------
  # Reactive table 1: PROGRESS  ----------------
  output$saved_annotations <- shiny::renderTable({
    print("- Rendering table 1")
    
    table_output <- reactive_values$selected_cell_tags %>% 
      bind_rows(.id = "ucid_t.frame")  #%>% mutate(ucid = as.numeric(ucid_t.frame))
    
    if(nrow(table_output) > 0){
      table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
        mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
        left_join(select(d, ucid, pos, cellID))
    } else {
      table_output <- data.frame(message = "No annotations yet...")
    }
    
    table_output
  })
  
  # Reactive text 1: current UCID and t.frame  ----------------
  output$cell_ith <- shiny::renderText({
    ith_ucid <- d$ucid[reactive_values$ith_cell]
    
    tag_line <- paste(
      shiny::isolate(reactive_values$ith_cell),
      ith_ucid,
      # paste(unlist(reactive_values$selected_cell_tags[[ith_ucid]]), collapse = "; "),
      sep = ", "
    )
    
    # write(tag_line, file=tmp_csv_output,append=TRUE)
    
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
      print(paste("--", magick.cell$ucids))
    } else {
      # Output white if selection is empty
      print("-- Selection is empty")
      tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
    }
    
    tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
    
    list(src = tmpfile, 
         contentType = "image/jpeg")
  }, deleteFile=TRUE)
  
  # Reactive plot 1: user plot  ----------------
  output$plot <- shiny::renderPlot({
    print("- Rendering plot 1")
    
    ith_ucid <- as.character(d$ucid[reactive_values$ith_cell])
    ith_t.frame <- as.character(d$t.frame[reactive_values$ith_cell])
    
    print(paste("--", ith_ucid))
    
    ucid_data <- filter(cdata, ucid == ith_ucid)
    
    if(!is.null(tag_ggplot)){
      # Add data
      tag_ggplot_render <- tag_ggplot %+% ucid_data
      
      # Add current t.frame
      tag_ggplot_render <- tag_ggplot_render + geom_vline(xintercept = as.numeric(ith_t.frame), 
                                                          color = "black")
      
      # Add annotations
      table_output <- reactive_values$selected_cell_tags %>% 
        bind_rows(.id = "ucid_t.frame")  #%>% mutate(ucid = as.numeric(ucid_t.frame))
      
      if(nrow(table_output) > 0){
        table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
          mutate(ucid = as.integer(ucid), 
                 t.frame = as.integer(t.frame)) %>% 
          left_join(d[, c("ucid", "pos", "cellID", reactive_values$click_vars$yvar)]) %>% 
          filter(ucid == ith_ucid)
        
        table_output_longer <- table_output %>%
          select(-ucid, -cellID, -pos) %>%
          mutate_at(
            vars(one_of(names(cell_tags))),
            as.character
          ) %>%
          pivot_longer(-t.frame,
                       names_to = "categoria",
                       values_to = "valor",
                       values_drop_na = TRUE)
        
        tag_ggplot_render <- tag_ggplot_render + geom_vline(data = table_output_longer,
                                                            aes(xintercept = t.frame,
                                                                color = categoria,
                                                                text = valor
                                                            ), size = 2, linetype = 2) +
          theme(legend.position = "none")
        
        # table_output_longer_ith_ucid_yvals <- left_join(table_output_longer_ith_ucid,
        #                                        ucid_data[,c("t.frame", reactive_values$click_vars$yvar)]) 
          
        # print("-- Adding annotations to plot")
        # tag_ggplot_render <- tag_ggplot_render + 
        #   ggrepel::geom_label_repel(data = table_output_longer_ith_ucid_yvals,
        #                             aes(label = paste(categoria, valor, collapse = ": ")))
      }
      
      # Render
      print("-- Rendering plot")
      tag_ggplot_render
    }
  })
}


