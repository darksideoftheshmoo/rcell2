#' Shiny app server function for tagCell
#' Define server logic required to draw a histogram
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import shiny shinyjs formattable dplyr tidyr hexbin magick keys
#' @importFrom graphics polygon
tagCellServer <- function(input, output, session) {
  if(runtime_messages) print("tagCellServer 1: start")
  
  p <- paths
  d <- cdata %>% 
    dplyr::arrange(ucid, t.frame) %>% 
    {
      if(randomize_ucids){
        sample(x = split(., .$ucid), 
               size = length(unique(cdata$ucid))) %>% 
          bind_rows()
      } else {.}
    } %>% 
    mutate(ucid_t.frame = paste(ucid, t.frame, sep = "_")) %>% 
    mutate(
      cellID = as.integer(cellID),
      ucid = as.integer(ucid),
      t.frame = as.integer(t.frame)
    )
  
  
  ucid.unique <- unique(d$ucid)
  ucid.viewed <- setNames(rep(F, length(ucid.unique)), ucid.unique)
  
  if(!all.unique(d$ucid_t.frame)) 
    stop("tagCellServer error: the ucid-t.frame combination is not a primary key! (i.e it is not unique, check your data).")
  
  reactive_values <- shiny::reactiveValues(ith_cell = 1, # row index used to order dataframe rows
                                           ith_ucid = numeric(),
                                           i_line = 1,
                                           other_reactive_values = c(),
                                           selected_cell_tags = list(),
                                           click_vars = list(),
                                           ucid.viewed=ucid.viewed,
                                           next_key=FALSE,
                                           prev_key=FALSE,
                                           skip_key=FALSE,
                                           unskip_key=FALSE
                                           )
  
  # Restores selected tags list
  if(!is.null(tags.df)){
    previous.tags.list <- tags.df %>% 
      # remove entries with no t.frame info (i.e. viewed but not tagged cells)
      filter(!is.na(t.frame)) %>% 
      # reform the unique id colum
      unite("ucid_t.frame", ucid, t.frame, sep = "_") %>% 
      # cleanup columns, select only those on the cell_tags list
      .[c("ucid_t.frame", names(cell_tags))] %>% 
      # reform the tags list
      split(~ucid_t.frame) %>% 
      # reform list structure, removing NA entries
      lapply(function(.tags){
        .tags <- as.list(.tags[1,-1])
        .tags <- .tags[!sapply(.tags, is.na)]
        return(.tags)
      }) %>% 
      # ensure no empty list items
      {.[lapply(., length) > 0]}
  } else {
    previous.tags.list <- NULL
  }
  # Restore previous tagging reactive list
  if(!is.null(previous.tags.list)) reactive_values$selected_cell_tags <- previous.tags.list
  
  ### KEYS OBSERVER   ----------------
  observeEvent(input$keys, {
    if(runtime_messages) print("tagCellServer 12: keys observer fired")
    
    # reactive_values$pressed_key <- input$keys
    pressed_key <- input$keys
    
    if(runtime_messages) print(paste("tagCellServer 12: Pressed keys:", pressed_key, collapse = " "))
    
    switch (pressed_key,
            `left` = {
              reactive_values$prev_key <- !reactive_values$prev_key
            },
            `right` = {
              reactive_values$next_key <- !reactive_values$next_key
            },
            `shift+left` = {
              reactive_values$unskip_key <- !reactive_values$unskip_key
            },
            `shift+right` = {
              reactive_values$skip_key <- !reactive_values$skip_key
            }
    )
  })
  
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
    if(runtime_messages) print("tagCellServer 2: renderUI update")
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
    if(runtime_messages) print("tagCellServer 3: plot_click event observer")
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
  toListen.next <- reactive({
    list(input$next_cell,
         reactive_values$next_key)
  })
  shiny::observeEvent(
    ignoreInit = T,
    eventExpr = toListen.next(),
    handlerExpr = {
      if(runtime_messages) print("tagCellServer 4: next_cell event observer")
      if(debug_messages) print("- Next cell requested, saving current tags...")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_cell")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid for that cell
      
      # Mark current ucid as viewed
      ith_ucid2 <- as.character(d$ucid[reactive_values$ith_cell])  # Get ucid for that cell
      reactive_values$ucid.viewed[names(ucid.viewed) == ith_ucid2] <- TRUE
      
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
  toListen.prev <- reactive({
    list(input$prev_cell,
         reactive_values$prev_key)
  })
  shiny::observeEvent(
    ignoreInit = T,
    # eventExpr = input$prev_cell,
    eventExpr = toListen.prev(),
    handlerExpr = {
      if(runtime_messages) print("tagCellServer 5: prev_cell event observer")
      if(debug_messages) print("- Previous cell requested...")
      shinyjs::disable("prev_cell")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                               # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell]) # Get ucid for that cell
      
      # Mark current ucid as viewed
      ith_ucid2 <- as.character(d$ucid[reactive_values$ith_cell])  # Get ucid for that cell
      reactive_values$ucid.viewed[names(ucid.viewed) == ith_ucid2] <- TRUE
      
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
  toListen.skip <- reactive({
    list(input$next_ucid,
         reactive_values$skip_key)
  })
  shiny::observeEvent(
    ignoreInit = T,
    # eventExpr = input$next_ucid,
    eventExpr = toListen.skip(),
    handlerExpr = {
      if(runtime_messages) print("tagCellServer 6: next_ucid event observer")
      if(debug_messages) print("- Next ucid requested, saving current tags...")
      shinyjs::disable("next_ucid")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_cell")
      shinyjs::disable("prev_cell")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid_t.frame for that cell
      
      # Mark current ucid as viewed
      ith_ucid2 <- as.character(d$ucid[reactive_values$ith_cell])  # Get ucid for that cell
      reactive_values$ucid.viewed[names(ucid.viewed) == ith_ucid2] <- TRUE
      
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
          ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags  # Save the tag list to a UCID name element in a reactive values list.
      
      # Skip to the next UCID
      ucid.oi <- as.character(d$ucid[reactive_values$ith_cell])   # Get bare ucid for the current cell
      ucid.oi.index <- match(ucid.oi, ucid.unique)                # Get it's index in the unique ucids array
      ucid.next <- ucid.unique[ucid.oi.index + 1]                 # Get the next ucid in that array
      ucid.next.index <- match(ucid.next, d$ucid)                 # And get the next ucid's row index
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
  toListen.unskip <- reactive({
    list(input$prev_ucid,
         reactive_values$unskip_key)
  })
  shiny::observeEvent(
    ignoreInit = T,
    # eventExpr = input$prev_ucid,
    eventExpr = toListen.unskip(),
    handlerExpr = {
      if(runtime_messages) print("tagCellServer 7: prev_ucid event observer")
      if(debug_messages) print("- Previous ucid requested...")
      shinyjs::disable("prev_ucid")
      shinyjs::disable("next_ucid")
      shinyjs::disable("prev_cell")
      shinyjs::disable("next_cell")
      shinyjs::disable("moreControls")
      
      ith_cell <- reactive_values$ith_cell                               # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell]) # Get ucid for that cell
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
      
      # Mark current ucid as viewed
      ith_ucid2 <- as.character(d$ucid[reactive_values$ith_cell])  # Get ucid for that cell
      reactive_values$ucid.viewed[names(ucid.viewed) == ith_ucid2] <- TRUE
      
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
    if(runtime_messages) print("tagCellServer 8: ith_cell reactive value observer")
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
  
  # SESSION END OBSERVER: EXIT  ----------------
  # This code will be run after the client has disconnected
  # https://stackoverflow.com/a/33236190/11524079
  # You'll have to use isolate to access reactiveValues (such as session) in a non reactive context
  # https://stackoverflow.com/a/49384590/11524079
  session$onSessionEnded(function() isolate({
    writeLines("\n- Quit event fired")
    
    if(debug_messages) print(paste("-- Saving progress to file:", tmp_output_file))
    
    # Save current annotations
    ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
    ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid_t.frame for that cell
    
    if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
    ith_cell_tags <- list()
    for(tag_group in 1:length(names(cell_tags)))      # For each tag group
      input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
      ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
    # Save the tag list to a UCID name element in a reactive values list.
    reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags
    
    # Mark current ucid as viewed
    ith_ucid2 <- as.character(d$ucid[reactive_values$ith_cell])  # Get ucid for that cell
    reactive_values$ucid.viewed[names(ucid.viewed) == ith_ucid2] <- TRUE
    # Generate viewed ucids DF
    viewed_ucids <- data.frame(ucid = as.integer(names(reactive_values$ucid.viewed)),
                               viewed = unname(reactive_values$ucid.viewed))
    
    # Bind current tags
    table_output <- reactive_values$selected_cell_tags %>% 
      bind_rows(.id = "ucid_t.frame")
    
    if(nrow(table_output) > 0){
      
      table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame"), convert = T) %>% 
        mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
        right_join(viewed_ucids)
      
    } else {
      table_output <- viewed_ucids
    }
    
    readr::write_csv(table_output, path = tmp_output_file)
    
    if(debug_messages) print("-- Returning progress to output:")
    stopApp(table_output)
  }))
  
  # BUTTON OBSERVER 3: EXIT  ----------------
  observeEvent(
    # Acción al apretar el botón de cerrar la app
    eventExpr = input$quit,
    handlerExpr = {
      writeLines("\n- Quit event fired")

      if(debug_messages) print(paste("-- Saving progress to file:", tmp_output_file))
      
      # Save current annotations
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid_t.frame for that cell
      
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
        ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      # Save the tag list to a UCID name element in a reactive values list.
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags
      
      # Mark current ucid as viewed
      ith_ucid2 <- as.character(d$ucid[reactive_values$ith_cell])  # Get ucid for that cell
      reactive_values$ucid.viewed[names(ucid.viewed) == ith_ucid2] <- TRUE
      # Generate viewed ucids DF
      viewed_ucids <- data.frame(ucid = as.integer(names(reactive_values$ucid.viewed)),
                                 viewed = unname(reactive_values$ucid.viewed))
      
      # Bind current tags
      table_output <- reactive_values$selected_cell_tags %>% 
        bind_rows(.id = "ucid_t.frame")
      
      if(nrow(table_output) > 0){
        
        table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame"), convert = T) %>% 
          mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
          right_join(viewed_ucids)
        
      } else {
        table_output <- viewed_ucids
      }
      
      readr::write_csv(table_output, path = tmp_output_file)
      
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
      
      # Save current annotations
      ith_cell <- reactive_values$ith_cell                                 # Get the current reactive cell number
      ith_ucid <- as.character(d$ucid_t.frame[reactive_values$ith_cell])   # Get ucid_t.frame for that cell
      
      if(debug_messages) print(paste("-- Saving tag selection for current cell with row index:", ith_cell))
      ith_cell_tags <- list()
      for(tag_group in 1:length(names(cell_tags)))      # For each tag group
        input[[names(cell_tags)[tag_group]]] ->         # Get the currently selected values array
        ith_cell_tags[[names(cell_tags)[tag_group]]]  # Store it in a list element appropriately named 
      # Save the tag list to a UCID name element in a reactive values list.
      reactive_values$selected_cell_tags[[ith_ucid]] <- ith_cell_tags
      
      # Mark current ucid as viewed
      ith_ucid2 <- as.character(d$ucid[reactive_values$ith_cell])  # Get ucid for that cell
      reactive_values$ucid.viewed[names(ucid.viewed) == ith_ucid2] <- TRUE
      # Generate viewed ucids DF
      viewed_ucids <- data.frame(ucid = as.integer(names(reactive_values$ucid.viewed)),
                                 viewed = unname(reactive_values$ucid.viewed))
      
      # Prepare output
      table_output <- reactive_values$selected_cell_tags %>% 
        bind_rows(.id = "ucid_t.frame")  #%>% mutate(ucid = as.numeric(ucid_t.frame))
      
      if(nrow(table_output) > 0){
        table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
          mutate(ucid = as.integer(ucid), t.frame = as.integer(t.frame)) %>% 
          right_join(viewed_ucids)
        
        showNotification(paste("-- Saving progress to file:", tmp_output_file), duration = 4, type = "message")
        
      } else {
        table_output <- viewed_ucids
        showNotification(paste("-- No annotations yet, nothing was saved."), duration = 4, type = "message")
      }
      
        # Write output
        readr::write_csv(table_output, path = tmp_output_file, append = F)
    }
  )
  
  ### OUTPUT OBSERVERS and RENDERERS  ----------------
  # Reactive TABLE 1: annotation progress  ----------------
  ## Input reactive_values: $selected_cell_tags
  ## Output reactive_values: 
  ## Isolated reactive_values: $selected_cell_tags
  output$saved_annotations <- shiny::renderTable({
    if(runtime_messages) print("tagCellServer 9: selected_cell_tags reactive table observer")
    if(debug_messages) print("- Rendering table 1")
    
    table_output <- reactive_values$selected_cell_tags %>% 
      data.table::rbindlist(idcol = "ucid_t.frame", fill=TRUE)
      # bind_rows(.id = "ucid_t.frame")  #%>% mutate(ucid = as.numeric(ucid_t.frame))
    
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
  
  # Reactive IMAGE 2: cell strips  ----------------
  output$pics2 <- shiny::renderImage({
    if(runtime_messages) print("tagCellServer 10: ith_cell reactive image observer 2")
    if(debug_messages) print("- Rendering image 2")
    
    # Make the image match the plot's width
    if(is.null(cell_resize)){
      output_plot_width <- session$clientData$output_plot_width
      cell_resize <- output_plot_width/length(tag_channels_select)
    }
    # session$clientData$output_plot_height
    
    if(nrow(d) > 0 & max.frames > 0) {
      if(debug_messages) print("-- Selection not empty: magick strip!")
      
      # d <- cdata.continuous
      # reactive_values <- list(ith_cell = 2)
      # plot_width <- 300
      # plot_height <- 200
      # tag_box_size = 50
      # equalize_images = F
      # normalize_images = F
      # tag_channels_select=c("BF", "BF.out")
      
      ucid.selected <- d[reactive_values$ith_cell,"ucid",drop=T]
      frame.selected <- d[reactive_values$ith_cell,"t.frame",drop=T]
      
      avail.frames <- filter(d, ucid == ucid.selected) %>% with(t.frame) %>% sort()
      
      # Reorder available frames by distance to current frame
      closest.frames <- avail.frames[order(abs(avail.frames - frame.selected))]
      # Truncate the array to max.frames length (or total length)
      closest.frames <- closest.frames[1:(min(length(closest.frames), max.frames))]
      # Sort the frames, just 'cause:
      frame.range <- sort(closest.frames)
      
      cell.strip <- 
        cellStrips(cdata = cdata %>% filter(ucid == ucid.selected,
                                            t.frame %in% frame.range),
                   paths = paths,
                   # cell_resize=cell_resize,
                   n.cells = max.frames,
                   ch=tag_channels_select, 
                   equalize_images = equalize_images,
                   normalize_images = normalize_images,
                   boxSize = tag_box_size, 
                   ) %>% 
        # Unlist the output
        .[[1]]
      
      
      plot_width <- session$clientData$output_plot_width
      plot_height <- session$clientData$output_plot_height
      
      tmpimage <- cell.strip %>% 
        magick::image_resize(geometry = geometry_size_pixels(width = plot_width))
      
      if(debug_messages) print(paste("--", magick.cell$ucids))
    } else {
      # Output white if selection is empty
      if(debug_messages) print("-- Selection is empty")
      tmpimage <- magick::image_blank(150,10,color = "white") %>% image_annotate(text = "(cell strip placeholder)")
    }
    
    tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
    
    list(src = tmpfile, 
         contentType = "image/jpeg")
  }, deleteFile=TRUE)
  
  # Reactive IMAGE 1: cell magick  ----------------
  output$pics <- shiny::renderImage({
    if(runtime_messages) print("tagCellServer 10: ith_cell reactive image observer")
    if(debug_messages) print("- Rendering image 1")
    
    # Make the image match the plot's width
    if(is.null(cell_resize)){
      output_plot_width <- session$clientData$output_plot_width
      cell_resize <- output_plot_width/length(tag_channels_select)
    }
    # session$clientData$output_plot_height
    
    if(nrow(d) > 0) {
      if(debug_messages) print("-- Selection not empty: magick!")
      # cdata.selected <- d[d$ucid == d$ucid[reactive_values$ith_cell],]
      cdata.selected <- d[reactive_values$ith_cell,]
      magick.cell <-  magickCell(cdata = cdata.selected, 
                                 p,
                                 # ch=input$image_channel, 
                                 # cell_resize=cell_resize,
                                 cell_resize=cell_resize,
                                 ch=tag_channels_select, 
                                 equalize_images = equalize_images,
                                 normalize_images = normalize_images,
                                 boxSize = tag_box_size, 
                                 return_single_imgs = T, 
                                 return_ucid_df = T)
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
  
  # Reactive PLOT: user plot  ----------------
  output$plot <- shiny::renderPlot({
    if(runtime_messages) print("tagCellServer 11: ith_cell reactive user plot observer")
    if(debug_messages) print("- Rendering plot 1")
    
    ith_ucid <- as.character(d$ucid[reactive_values$ith_cell])
    ith_t.frame <- as.character(d$t.frame[reactive_values$ith_cell])
    
    if(debug_messages) print(paste("-- Current reactive_values$ith_cell:", reactive_values$ith_cell))
    if(debug_messages) print(paste("-- Current ith_ucid:", ith_ucid))
    
    ucid_data <- filter(cdata, ucid %in% ith_ucid)
    
    if(is.null(tag_ggplot)){
      if(runtime_messages) print("tagCellServer 11: no tag_ggplot object provided, rendering defult plot!")
      tag_ggplot <- ggplot() +
        geom_line(aes(x=t.frame, y=a.tot)) + 
        geom_vline(xintercept = as.numeric(ith_t.frame), color = "black", linetype=2)
      tag_ggplot_render <- tag_ggplot %+% ucid_data
    } else {
      # Add data
      tag_ggplot_render <- tag_ggplot %+% ucid_data
      
      # Add current t.frame
      if(debug_messages) print(paste("-- adding geom_vline to plot:"))
      tag_ggplot_render <- tag_ggplot_render + geom_vline(xintercept = as.numeric(ith_t.frame), 
                                                          color = "black", linetype=2)
      
      # Add annotations
      if(debug_messages) print(paste("-- generating selected_cell_tags"))
      # if(debug_messages) print(names(reactive_values$selected_cell_tags))
      # if(debug_messages) print(data.table::rbindlist(reactive_values$selected_cell_tags, idcol = "ucid_t.frame"))
      table_output <- reactive_values$selected_cell_tags %>% 
        data.table::rbindlist(idcol = "ucid_t.frame", fill=TRUE)
        # bind_rows(.id = "ucid_t.frame")  #%>% mutate(ucid = as.numeric(ucid_t.frame))
      
      # if(debug_messages) print(table_output)
      if(nrow(table_output) > 0){
        if(debug_messages) print(paste("-- found annotations for plot, processing..."))
        table_output <- separate(table_output, ucid_t.frame, c("ucid", "t.frame")) %>% 
          mutate(ucid = as.integer(ucid), 
                 t.frame = as.integer(t.frame)) %>% 
          left_join(d[, c("ucid", "pos", "cellID", reactive_values$click_vars$yvar)]) %>% 
          filter(ucid == ith_ucid)
        
        table_output_longer <- table_output %>%
          select(-ucid, -cellID, -pos) %>%
          mutate_at(
            dplyr::vars(tidyselect::any_of(names(cell_tags))),
            as.character
          ) %>%
          pivot_longer(-t.frame,
                       names_to = "categoria",
                       values_to = "valor",
                       values_drop_na = TRUE)
        
        if(debug_messages) print(paste("-- adding geom_vline to plot:"))
        tag_ggplot_render <- tag_ggplot_render + geom_vline(data = table_output_longer,
                                                            aes(xintercept = t.frame,
                                                                color = interaction(categoria, valor, sep = ": "),
                                                                # linetype = valor,
                                                                # text = paste(categoria, valor, sep = ": ")
                                                                ),
                                                            size = 2, linetype = 2) +
          # ggplot2::guides(text = FALSE) +  # http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software
          theme(legend.position = "bottom", legend.title = element_blank())
      }
    }
    
    # Render
    if(debug_messages) print("-- Rendering plot")
    tag_ggplot_render
  })
}


