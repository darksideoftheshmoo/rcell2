#' Shiny app server function for shinyCell
#' Define server logic required to draw a histogram
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import shiny dplyr tidyr
plotAppServer <-
  function(input, output, session) {
    
    # Get plot variables ####
    .variables <- c(xvar = quo_name(user_plot$layers[[1]]$mapping$x),
                    yvar = quo_name(user_plot$layers[[1]]$mapping$y))
    
    
    if(debug_messages) print("-- Plot variables")
    if(debug_messages) print(.variables)
    
    output$scatterplot <- renderPlot({
      user_plot
    })
    
    # Table for plot brush ####
    output$selection_table <- shiny::renderDataTable({
      
      if(!is.null(input$scatterplot_brush$xmin)){
        
        brush_limits <- c(xmin = input$scatterplot_brush$xmin,
                          xmax = input$scatterplot_brush$xmax,
                          ymin = input$scatterplot_brush$ymin,
                          ymax = input$scatterplot_brush$ymax)
        
        if(debug_messages) print("-- Brush limits")
        if(debug_messages) print(brush_limits)
        
        
        # p.data <- ggplot2::layer_data(p, 1)
        # plot.data <- user_plot$data
        
        # selector <- 
        #   as.numeric(user_data[, .variables["xvar"]]) >= brush_limits["xmin"] &
        #   as.numeric(user_data[, .variables["xvar"]]) <= brush_limits["xmax"] &
        #   as.numeric(user_data[, .variables["yvar"]]) <= brush_limits["ymax"] &
        #   as.numeric(user_data[, .variables["yvar"]]) >= brush_limits["ymin"]
        # # user_data[selector, ]
        
        selected <- filter(user_plot$data,
                           as.numeric( !!sym(.variables[["xvar"]]) ) >= brush_limits[["xmin"]] &
                           as.numeric( !!sym(.variables[["xvar"]]) ) <= brush_limits[["xmax"]] &
                           as.numeric( !!sym(.variables[["yvar"]]) ) <= brush_limits[["ymax"]] &
                           as.numeric( !!sym(.variables[["yvar"]]) ) >= brush_limits[["ymin"]])
      } else {
        selected <- data.frame()
      }
      # Return selected points table
      selected
    }, 
    options = list(scrollX = TRUE, pageLength = 10))
    
    # Quit event ####
    observeEvent(
      eventExpr = input$quit,
      handlerExpr = {
        
        if(!is.null(input$scatterplot_brush$xmin)){
          
          brush_limits <- c(xmin = input$scatterplot_brush$xmin,
                            xmax = input$scatterplot_brush$xmax,
                            ymin = input$scatterplot_brush$ymin,
                            ymax = input$scatterplot_brush$ymax)
          
          selected <- filter(user_plot$data,
                             as.numeric( !!sym(.variables[["xvar"]]) ) >= brush_limits[["xmin"]] &
                               as.numeric( !!sym(.variables[["xvar"]]) ) <= brush_limits[["xmax"]] &
                               as.numeric( !!sym(.variables[["yvar"]]) ) <= brush_limits[["ymax"]] &
                               as.numeric( !!sym(.variables[["yvar"]]) ) >= brush_limits[["ymin"]])
          
        } else {
          selected <- data.frame()
        }
        
        stopApp(selected)
      }
    )
  } # end server definition
