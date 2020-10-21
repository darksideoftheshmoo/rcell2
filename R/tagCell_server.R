#' Shiny app server function for tagCell
#' Define server logic required to draw a histogram
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import shiny shinyjs formattable dplyr tidyr hexbin magick
#' @importFrom graphics polygon
tagCellServer <- function(input, output, session) {
  
  if(is.null(tmp_output_file)){
    tmp_output_file <- tempfile(tmpdir = "./", fileext = ".txt")
  }
  if(debug_messages) print(paste("Appending tags to tempfile:", tmp_output_file)) # find: ^([\s\t]+)print replace: \1if(debug_messages) print
  dir.create(dirname(normalizePath(tmp_output_file)), recursive = T)
  # write("", file=tmp_output_file)
  
  d <- cdata %>% 
    dplyr::arrange(ucid, t.frame) %>% 
    {if(randomize_ucids) 
      sample(x = split(., .$ucid), 
             size = length(unique(cdata$ucid))) %>% 
        bind_rows() else .} %>% 
    mutate(ucid_t.frame = paste(ucid, t.frame, sep = "_")) %>% 
    mutate(
      cellID = as.integer(cellID),
      ucid = as.integer(ucid),
      t.frame = as.integer(t.frame)
    )
  p <- paths
  
  ucid.unique <- unique(d$ucid)
  
  reactive_values <- shiny::reactiveValues(ith_cell = 1, # row index used to order dataframe rows
                                           ith_ucid = numeric(),
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
    if(debug_messages) print("- Plot clicked")
    shinyjs::disable("next_cell")
    shinyjs::disable("prev_cell")
    shinyjs::disable("prev_ucid")
    shinyjs::disable("next_ucid")
    shinyjs::disable("moreControls")
    
    closest_to <- function(closest_to, from_array){
      if(debug_messages) print(paste(">> Closest to:", closest_to, "from array:", paste(from_array, collapse = ", ")))
      index <- which.min(abs(from_array - closest_to))
      return(from_array[index])
    }
    
    if(debug_messages) print(paste("-- Clicked point:", input$plot_click$x))
    
    ith_cell <- reactive_values$ith_cell                                        # Get the current reactive cell number
    ith_ucid_t.frame <- as.character(d$ucid_t.frame[reactive_values$ith_cell])  # Get ucid_t.frame for that cell
    ith_ucid <- as.character(d$ucid[reactive_values$ith_cell])                  # Get ucid for that cell
    ith_t.frame <- as.character(d$t.frame[reactive_values$ith_cell])            # Get t.frame for that cell
    
    click_t.frame <- closest_to(from_array = d[d$ucid == ith_ucid,]$t.frame,
                                closest_to = input$plot_click$x)
    
    if(debug_messages) print(paste("-- Clicked t.frame:", click_t.frame))
      
    if(click_t.frame == ith_t.frame){
      if(debug_messages) print("-- t.frame unchanged")
      shinyjs::enable("next_cell")
      shinyjs::enable("prev_cell")
      shinyjs::enable("prev_ucid")
      shinyjs::enable("next_ucid")
      shinyjs::enable("moreControls")
    } else {
      if(debug_messages) print("-- t.frame changed, saving tag selection")
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
        ith_cell_tags[[names(cell_tags)[tag_group]]]    # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid_t.frame]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      next_ith_cell <- which(d$t.frame == click_t.frame & d$ucid == ith_ucid)  # Get the index row for the clicked t.frame
      reactive_values$ith_cell <- next_ith_cell                                # Update the ith_cell reactive value
    }
  })
  
  # BUTTON 1.1: NEXT  ----------------
  ## Input: input$next_cell
  ## Output reactive_values: $selected_cell_tags $ith_cell
  shiny::observeEvent(
    eventExpr = input$next_cell,
    handlerExpr = {
      if(debug_messages) print("- Next cell requested, saving current tags...")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_cell")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid for that cell
      
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
          ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      # reactive_values$ith_ucid <- as.character(d$ucid[reactive_values$ith_cell])  # Save the ucid
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags             # Save the tag list to a UCID name element in a reactive values list.
      
      if(debug_messages) print(paste("-- Next cell row index:", ith_cell + 1))
      
      # Handle previous > total
      if(ith_cell == nrow(cdata)){
        showNotification("There is no next cell, staying at the current one.", type = "warning")
        shinyjs::delay(300, expr = {
          shinyjs::enable("prev_cell")
          shinyjs::enable("next_cell")
          shinyjs::enable("prev_ucid")
          shinyjs::enable("next_ucid")
          shinyjs::enable("moreControls")
        })
      } else {
        reactive_values$ith_cell <- ith_cell + 1                     # Update the ith_cell reactive value
      }
    })
  
  # BUTTON 1.2: PREVIOUS  ----------------
  shiny::observeEvent(
    eventExpr = input$prev_cell,
    handlerExpr = {
      if(debug_messages) print("- Previous cell requested...")
      shinyjs::disable("prev_cell")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                               # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell]) # Get ucid for that cell
      
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
          ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      if(debug_messages) print(paste("-- Next cell row index:", ith_cell - 1))
      
      # Handle previous < 1
      if(ith_cell == 1){
        if(debug_messages) print("-- There is no previous cell, staying at the current one.")
        showNotification("There is no previous cell, staying at the current one.", type = "warning")
        shinyjs::delay(300, expr = {
          shinyjs::enable("prev_cell")
          shinyjs::enable("next_cell")
          shinyjs::enable("prev_ucid")
          shinyjs::enable("next_ucid")
          shinyjs::enable("moreControls")
        })
      } else {
        if(debug_messages) print(paste("-- Moving on to cell with row index:", ith_cell - 1))
        reactive_values$ith_cell <- ith_cell - 1                   # Update the ith_cell reactive value
      }
    })
  
  # BUTTON 2.1: NEXT UCID  ----------------
  shiny::observeEvent(
    eventExpr = input$next_ucid,
    handlerExpr = {
      if(debug_messages) print("- Next ucid requested, saving current tags...")
      shinyjs::disable("next_ucid")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_cell")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid_t.frame for that cell
      
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
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
      if(debug_messages) print(paste("-- Next cell row index:", ucid.next.index))
      
      # Handle next > total
      if(ucid.oi.index >= length(ucid.unique)){
        showNotification("--- There is no next ucid, staying at the current one.", type = "warning")
        shinyjs::delay(300, expr = {
          shinyjs::enable("prev_cell")
          shinyjs::enable("next_cell")
          shinyjs::enable("prev_ucid")
          shinyjs::enable("next_ucid")
          shinyjs::enable("moreControls")
        })
      } else {
        showNotification("--- Moving to next ucid.", duration = 1)
        reactive_values$ith_cell <- ucid.next.index               # Update the ith_cell reactive value
      }
    })
  
  # BUTTON 2.2: PREVIOUS UCID ----------------
  shiny::observeEvent(
    eventExpr = input$prev_ucid,
    handlerExpr = {
      if(debug_messages) print("- Previous ucid requested...")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("prev_cell")
      shinyjs::disable("next_cell")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                               # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell]) # Get ucid for that cell
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
      
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
        ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      # Skip to the previous UCID
      ucid.oi <- as.character(d$ucid[reactive_values$ith_cell])   # Get bare ucid for the current cell
      ucid.oi.index <- match(ucid.oi, ucid.unique)                # Get it's index in the distinct ucid list 
      ucid.next <- ucid.unique[ucid.oi.index - 1]                 # Get the previous unique ucid
      ucid.next.index <- match(ucid.next, d$ucid)                 # Get the previous ucid's row index
      if(debug_messages) print(paste("-- Next cell row index:", ucid.next.index))
      
      # Handle next > total
      if(ucid.oi.index == 1){
        showNotification("--- There is no previous ucid, staying at the current one.", type = "warning")
        shinyjs::delay(300, expr = {
          shinyjs::enable("prev_cell")
          shinyjs::enable("next_cell")
          shinyjs::enable("prev_ucid")
          shinyjs::enable("next_ucid")
          shinyjs::enable("moreControls")
        })
      } else {
        showNotification("--- Moving to previous ucid.", duration = 1)
        reactive_values$ith_cell <- ucid.next.index               # Update the ith_cell reactive value
      }
    })
  
  # SIDE EFFECTS FOR PREV/NEXT BUTTON  ----------------
  ## Input reactive_values: $ith_cell
  ## Output reactive_values: 
  ## Isolated reactive_values: $selected_cell_tags
  shiny::observe({
    shinyjs::disable("prev_cell")
    shinyjs::disable("next_cell")
    shinyjs::disable("prev_ucid")
    shinyjs::disable("next_ucid")
    shinyjs::disable("moreControls")
    
    ith_cell <- reactive_values$ith_cell
    ith_ucid <- as.character(d$ucid_t.frame[ith_cell])
    selected_cell_tags <- isolate(reactive_values$selected_cell_tags)
    if(debug_messages) print(paste("-- Updating tag selection for next or previous cell with row index:", ith_cell))
    
    if(ith_ucid %in% names(selected_cell_tags)){
      if(debug_messages) print("--- UCID tag found")
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
      if(debug_messages) print("--- UCID not yet tagged")
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

      if(debug_messages) print(paste("-- Saving progress to file:", tmp_output_file))
      
      table_output <- reactive_values$selected_cell_tags %>% bind_rows(.id = "ucid_t.frame")
      
      if(nrow(table_output) > 0){
        
        table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
          mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
          left_join(unique(select(d, ucid, pos, cellID)))
        
        table_output %>% readr::write_csv(path = tmp_output_file)
        
      } else {
        
        table_output <- data.frame()
        
      }
      
      if(debug_messages) print("-- Returning progress to output:")
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
  # Reactive TABLE 1: annotation progress  ----------------
  ## Input reactive_values: $selected_cell_tags
  ## Output reactive_values: 
  ## Isolated reactive_values: $selected_cell_tags
  output$saved_annotations <- shiny::renderTable({
    if(debug_messages) print("- Rendering table 1")
    
    table_output <- reactive_values$selected_cell_tags %>% 
      bind_rows(.id = "ucid_t.frame")  #%>% mutate(ucid = as.numeric(ucid_t.frame))
    
    if(nrow(table_output) > 0){
      table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
        mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
        left_join(unique(select(d, ucid, pos, cellID)))
    } else {
      table_output <- data.frame(message = "No annotations yet...")
    }
    
    table_output
  })
  
  # Reactive TEXT 1: current UCID and t.frame  ----------------
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
  
  # Reactive TEXT 2: hover info  ----------------
  output$hover_info <- renderText({
    tryCatch(
      expr = {
        paste0("x=", signif(input$plot_hover$x, 3), "\t y=", signif(input$plot_hover$y, 3))
      },
      error = function(cond) return("Waiting for valid mouse hover...")
    )
  })
  
  # Reactive IMAGE 1: magickCell  ----------------
  output$pics <- shiny::renderImage({
    if(debug_messages) print("- Rendering image 1")
    
    if(nrow(d) > 0) {
      if(debug_messages) print("-- Selection not empty: magick!")
      # cdata.selected <- d[d$ucid == d$ucid[reactive_values$ith_cell],]
      cdata.selected <- d[reactive_values$ith_cell,]
      magick.cell <-  magickCell(cdata = cdata.selected, 
                                 p,
                                 # ch=input$image_channel, 
                                 cell_resize=cell_resize,
                                 ch=tag_channels_select, 
                                 n = n_max, 
                                 equalize_images = equalize_images,
                                 normalize_images = normalize_images,
                                 boxSize = tag_box_size, 
                                 return_single_imgs = T, ...)
      tmpimage <- magick.cell$img
      if(debug_messages) print(paste("--", magick.cell$ucids))
    } else {
      # Output white if selection is empty
      if(debug_messages) print("-- Selection is empty")
      tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
    }
    
    tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
    
    list(src = tmpfile, 
         contentType = "image/jpeg")
  }, deleteFile=TRUE)
  
  # Reactive PLOT 1: user plot  ----------------
  output$plot <- shiny::renderPlot({
    if(debug_messages) print("- Rendering plot 1")
    
    ith_ucid <- as.character(d$ucid[reactive_values$ith_cell])
    ith_t.frame <- as.character(d$t.frame[reactive_values$ith_cell])
    
    if(debug_messages) print(paste("--", ith_ucid))
    
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
                                                                color = interaction(categoria, valor, sep = ": "),
                                                                # linetype = valor,
                                                                text = paste(categoria, valor, sep = ": ")),
                                                            size = 2, linetype = 2) +
          # ggplot2::guides(text = FALSE) +  # http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software
          theme(legend.position = "bottom", legend.title = element_blank())
        
        # table_output_longer_ith_ucid_yvals <- left_join(table_output_longer_ith_ucid,
        #                                        ucid_data[,c("t.frame", reactive_values$click_vars$yvar)]) 
          
        # print("-- Adding annotations to plot")
        # tag_ggplot_render <- tag_ggplot_render + 
        #   ggrepel::geom_label_repel(data = table_output_longer_ith_ucid_yvals,
        #                             aes(label = paste(categoria, valor, collapse = ": ")))
      }
      
      # Render
      if(debug_messages) print("-- Rendering plot")
      tag_ggplot_render
    }
  })
}


