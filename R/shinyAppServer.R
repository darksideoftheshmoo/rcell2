#' Shiny app server function for shinyCell
#' Define server logic required to draw a histogram
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import shiny formattable dplyr tidyr hexbin magick
#' @importFrom graphics polygon
shinyAppServer <-
  function(input, output, session) {

    ### INITIALIZATION ###
    {
      # Reactive values definition
      values <- reactiveValues()
      
      # Initialize cdata
      cdata$filter <- T
      if(!{cdata$t.frame %>% unique() %>% length()} == 1){
        print("Initializing cell_unique_id_field to 'ucid_time'")
        cdata <- mutate(cdata, ucid_time = paste0(ucid, "_", t.frame))
        cdata.cell_unique_id_field <- "ucid_time"
      } else {
        print("Initializing cell_unique_id_field to 'ucid'")
        cdata.cell_unique_id_field <- "ucid"
      }
      values$cdata <- cdata
      
      # Initialize cfilter
      # values$cfilter <- data.frame(ucid = cdata[,"ucid"], filter = T)
      # if(length(filters) > 0) { values$stringFilters <- filters} else {values$stringFilters <- c()}  # Load initial filters if any
      # values$stringFiltersSelected <- c()

      # Initialize sample seed
      values$seed <- seed
      
      # Initialuze hover for nearby cells
      # values$points_nearby <- c()

      # Setup polygon
      pgnpts_empty <- data.frame(x = numeric(),
                                 y = numeric(),
                                 xvar = character(),
                                 yvar = character(),
                                 stringsAsFactors = F)

      rv <- reactiveValues(pgnpts = pgnpts_empty,
                           filters = filters,
                           brush_limits = c(),
                           others = c()
      )

      # Initialize output list, to fix: "  shinyCell : server: no visible binding for '<<-' assignment to ‘saved’"
      saved <- list()
      print("Initialization done!")
    }

    ### FILTER TAB OBSERVER ###
    {
    # Observador 01: cambios en filtros ----------------
      # Reads input: input$stringFilters input$truth_mode
      # Reads reactive: rv$filters isolate(values$cdata$filter)
      # Writes reactive: values$cdata values$cfilter
      observe({
        # Revampeado para usar bien geom_hex: https://stackanswers.net/questions/setting-hex-bins-in-ggplot2-to-same-size
        # Isolate from values$cdata, but not from input$stringFilters, rv$filters
        print("Filter update fired")

        # Si hay filtros tickeados, construir filtros y aplicarlos a cdata
        if(length(input$stringFilters) > 0) {
          print("-- Applying polygon filters selected")

          # Por ahora asumo que los filtros seleccionados tienen valores enteros en el nombre
          selectedPolygons <- as.numeric(input$stringFilters)

          # De la lista de filtros/dataframes, tomar solamente los seleccionados.
          selectedPolygons <- rv$filters[selectedPolygons]

          result <- polyFilterApply(polygon_df_list = selectedPolygons,
                                    cdata = cdata,
                                    truthMode = input$truth_mode, 
                                    cell_unique_id_field = cdata.cell_unique_id_field)  # important for time series

          # Update reactive cdata and cfilter
          # print(identical(values$cdata, result$cdata))  # Why??
          values$cdata <- result$cdata
          # values$cfilter <- result$cfilter

        } else {
          print("-- No polygon filters selected")

          # If no polygon filters are selected, reset that shit!
          values$cdata <- cdata
          # values$cfilter <- data.frame(ucid = cdata[,"ucid"], filter = T)
        }
        
        print(paste("-- Saving filter progress to:", filter_progress_file))
        saveRDS(object = rv$filters, file = filter_progress_file)
        print(paste("-- Remaining cells:", sum(isolate(values$cdata$filter))))  # Isolate from values$cdata

      }, priority = 11)
    }

    ### OUTPUT OBSERVERS ###
    {
      # Observador 02: Render plot event ----------------
      # Reads input: input$suspend_filters, input$ptype, input$x, input$y, input$stringFilters, input$overlay_polygons, input$facet input$plotDimY input$truth_mode
      # Reads reactive: isolate(values$cdata), rv$pgnpts, rv$filters
      # Writes reactive: -
      observe(priority = 1,
              label = "Plot observer",
              x = {
                output$scatterplot <- renderPlot({
                  print("Render plot event fired")

                  # Forcing update by truth criteria change
                  input$truth_mode

                  # Prepare stuff
                  plot.type <- input$ptype

                  positions <- rangeExpand(input$position, max(paths$pos))

                  if(input$suspend_filters){
                    d <- subset(isolate(values$cdata),  # isolated, plot is only updated by input$ change
                                pos %in% positions)
                  } else {
                    d <- subset(isolate(values$cdata),  # isolated, plot is only updated by input$ change
                                pos %in% positions & filter == TRUE)
                  }

                  print(paste0("-- Facet formula: ", input$facet))
                  facetVars <- getFacetVars(pdata, facetFormulaString = input$facet)
                  
                  print(paste("-- Rows left in filtered dataframe for plotting:", nrow(d)))
                  print(paste("-- Facet variables:", paste0(facetVars, collapse = ", ")))
                  print(paste("-- Facet vars in data:", paste0(facetVars %in% names(d), collapse = ", ")))

                  # Plot stuff
                  if(input$x == input$y){
                    plot.type <- "Density 1D"
                    print("X and Y are the same, forcing 1D density plot...")

                    p <- ggplot(d) + geom_density(aes(x = eval(parse(text=input$x))  #, y = ..level..
                                                      ),
                                                  position = "identity", alpha = .5) +
                      theme_minimal() + theme(legend.position="none") + xlab(input$x)

                    # TO-DO ####
                    # Density 1D does not play nice with brushes
                    # For these, the y-axis hould be meaningless and ignored
                    # Vertical brushing should be disabled
                    # Might be a good excuse to generalize drawing of polygons in plots where only one dimension is shared
                  }

                  if(plot.type == "Hex"){
                    print("-- Draw Hex plot...")

                    # Build HexPlot dataframe
                    # https://stackanswers.net/questions/setting-hex-bins-in-ggplot2-to-same-size
                    d <- hexPlotDf(d, facetVars = facetVars, varx = input$x, vary = input$y)
                    # Para fijar la relación de escala entre ejes en todos los facets
                    axisRatio = (max(d$x) - min(d$x)) / (max(d$y) - min(d$y))

                    p <- ggplot(data = d) +
                      geom_hex(aes(x = x, y = y, fill = counts), stat="identity") +
                      coord_equal(ratio = axisRatio) +
                      theme_bw()+ xlab(input$x) + ylab(input$y) +
                      scale_fill_continuous(type = "viridis", na.value = "#00000000")
                    #scale_fill_continuous(low = "grey80", high = "#000040", na.value = "#00000000")
                  }

                  if(plot.type == "Density"){
                    print("-- Draw density plot...")

                    # DensityPlot
                    p <- ggplot(data = d, aes(x = eval(parse(text=input$x)), y = eval(parse(text=input$y)))) +
                      geom_density2d() +
                      stat_density_2d(aes(fill = stat(level)), geom = "polygon", n = 2^7) +
                      theme(legend.position="none") + xlab(input$x) + ylab(input$y) +
                      scale_fill_continuous(type = "viridis", na.value = "#00000000")
                  }

                  if(plot.type == "Dots"){
                    print("-- Draw dot plot...")

                    # Dotplot
                    p <- ggplot(data = d, aes(x = eval(parse(text=input$x)), y = eval(parse(text=input$y)))) +
                      geom_point(alpha = 0.5) + {
                      if(input$x == "t.frame") 
                        geom_line(aes(group=ucid, color=factor(ucid)), alpha = 0.5) 
                      else 
                        geom_path(aes(group=ucid, color=factor(ucid)), alpha = 0.5)
                      } + 
                      geom_rug(col=grDevices::rgb(.5,0,0,alpha=.05)) +
                      xlab(input$x) + ylab(input$y)
                  }


                  # Overlay saved polygons, if any.
                  if(length(input$stringFilters) > 0 & input$overlay_polygons){
                    print("-- Overlaying saved polygon...")
                    pgnfilters.o <- bind_rows(rv$filters[as.numeric(input$stringFilters)],
                                            .id = "polygon")
                    # Get the relevant polygon filter points
                    pgnfilters <- subset(pgnfilters.o, xvar == input$x & yvar == input$y)
                    # Draw polygons even if axis have been swapped (i.e. a.tot VS fft.stat instead of fft.stat VS a.tot)
                    pgnfilters.swap.axis <- subset(pgnfilters.o, xvar == input$y & yvar == input$x) %>%
                      dplyr::rename(x = y,
                                    y = x)

                    p <- p + geom_polygon(data = bind_rows(pgnfilters, pgnfilters.swap.axis),
                                          aes(x = x, y = y,
                                              color = polygon,
                                              linetype = type),  # No se puede usar "fill", conflictua con "scale_fill" en Hex y Density
                                          size = 1,
                                          alpha = .1)

                    # TO DO
                    # Add labels to each polygon
                    # http://directlabels.r-forge.r-project.org/docs/lineplot/plots/chemqqmathscore.html
                    # https://stackoverflow.com/questions/29478152/plotting-text-labels-over-geom-polygon-data-in-ggmap-in-r
                  }

                  # Overlay newly drawn but not saved polygons
                  pgnpts <- subset(rv$pgnpts, xvar == input$x & yvar == input$y)
                  if(nrow(pgnpts) > 0){
                    p <- p +
                      geom_polygon(data = pgnpts,
                                   aes(x = x, y = y),
                                   color = "gray",
                                   alpha = .3) +
                      geom_point(data = pgnpts,
                                 aes(x = x, y = y),
                                 color = "red", size = 2)
                  }

                  # Draw facets if any
                  if(input$facet != "") {
                    facet <- parse(text=input$facet)

                    if(input$facet_grid) {
                      p <- p + facet_grid(eval(facet), scales = input$facet_scale)
                    } else {
                      p <- p + facet_wrap(eval(facet), scales = input$facet_scale)
                    }
                  }
                  
                  # Adjust plot limits
                  if(is.null(facets_scale_free)) {
                    if(plot.type != "Hex"){
                      if(plot.type == "Density 1D") {
                        p <- p + coord_cartesian(xlim = range(d[[input$x]]))
                      } else {
                        p <- p + coord_cartesian(xlim = range(d[[input$x]]), ylim = range(d[[input$y]]))
                      }
                    } else {
                      p <- p + coord_cartesian(xlim = range(d$x), ylim = range(d$y))
                    }
                  }

                  # Adjust theme
                  p <- p + 
                    theme_minimal() + 
                    theme(text = element_text(size=20),
                          legend.position = "none")

                  print("-- Plotting...")
                  p

                },
                height = input$plotDimY,
                width = input$plotDimX
                # https://shiny.rstudio.com/articles/client-data.html
                # width = "auto", height = "auto"
                )
              })
      # Observador 02.1: Session dimensions change ----------------
      # Reads input: input$plotDimX input$plotDimY session$clientData$output_scatterplot_width session$clientData$output_scatterplot_height
      # Reads reactive:
      # Writes reactive: input$plotDimX input$plotDimY
      observe(
        priority = 2,
        x = {
          # Nicer fix
          # https://shiny.rstudio.com/articles/client-data.html
          w1 <- isolate(input$plotDimX)
          w2 <- seq.int(500,2E3,50)[which.min(abs(session$clientData$output_scatterplot_width - seq.int(500,2E3,50)))]
          if(w1 != w2) updateSliderInput(session = session, inputId = "plotDimX", value = w2)

          h1 <- isolate(input$plotDimY)
          h2 <- seq.int(500,2E3,50)[which.min(abs(session$clientData$output_scatterplot_height - seq.int(500,2E3,50)))]
          if(h1 != h2) updateSliderInput(session = session, inputId = "plotDimY", value = h2 + 100)

          print("Session dimensions changed!")
          paste0(w1, " -w-> ", w2, " (", session$clientData$output_scatterplot_width, ")") %>% print()
          paste0(h1, " -h-> ", h2, " (", session$clientData$output_scatterplot_height, ")") %>% print()
        })

      # Observador 03: Brush observer ----------------
      # Reads input: input$scatterplot_brush$xmin, input$scatterplot_brush$xmax, input$scatterplot_brush$ymin, input$scatterplot_brush$ymax
      # Reads reactive: #values$cdata, isolate(rv$brush_limits)
      # Writes reactive: rv$brush_limits
      observe(priority = 5,
              label = "Brush observer",
              x = {
                print("Reactive values have been updated...")

                brush_limits <- c(input$scatterplot_brush$xmin,
                                  input$scatterplot_brush$xmax,
                                  input$scatterplot_brush$ymin,
                                  input$scatterplot_brush$ymax)

                # Redrawing the plot changes the last digit on the brush limits, this prevents futile firing of rv$brush_limits observers.
                if(!is.null(brush_limits)) brush_limits <- signif(brush_limits, 10)

                # uso isolate(rv$brush_limits) para no disparar el este mismo bloque al pedo (cuando se modifique rv$brush_limits)
                if(!identical(brush_limits, isolate(rv$brush_limits))){

                  #rv$brush_limits <- brush_limits
                  if( !is.null(input$scatterplot_brush$xmin) ) {
                    print("-- net change (not null)")
                    print(paste("---- Previous reactive brush: ", paste(isolate(rv$brush_limits), collapse = " ")))
                    print(paste("---- New reactive brush: ", paste(brush_limits, collapse = " ")))
                    rv$brush_limits <- brush_limits

                  } else print("-- brush change by null (not updating)")  # Si el nuevo brush está vacío, no hacer nada.

                } else print("-- no brush net change (not updating)")  # Si el nuevo brush es igual al anterior, no hacer nada.

                # Others
                # n <- nrow(values$cdata) # ??? lo comenté porque ni idea que hacia

              })


      # Reactive image 1 ----------------
      # Reads input: input$x, input$y, input$facet, input$facetinput$facet, input$facet_brush, input$scatterplot_brush, input$ch, input$position
      # Reads reactive: rv$brush_limits, values$cdata
      # Writes reactive: -
      output$pics <- renderImage({
        print("Rendering image 1: brush selection")

        brush_limits <- rv$brush_limits

        # Si hay algo seleccionado
        if(!is.null(brush_limits[1])){
          # Me quedo con las posiciones que me interesan
          positions <- rangeExpand(input$position, max(paths$pos))

          p <- subset(paths, pos %in% positions)
          d <- subset(values$cdata, pos %in% positions)
          if(!input$suspend_filters) d <- subset(d, filter == TRUE)

          # Y solamente con las filas que están en el brush
          d <- d[d[,input$x] >= brush_limits[1] &
                 d[,input$x] <= brush_limits[2] &
                 d[,input$y] >= brush_limits[3] &
                 d[,input$y] <= brush_limits[4] ,]

          # Si hay algo escrito en el facet formula field, filtrar el dataframe para mostrar solo las imagenes del facet seleccionado.
          if(input$facet != "" && input$facet_brush){
            print("-- Brush by facet mode ON")

            # Cada "eje" del facet tiene los posibles valores de una variable de pdata.

            # Primero conseguir las variables en la fórmula (debería estar primera la variable en el "eje y del facet")
            facetVars <- all.vars(eval(parse(text=input$facet)))
            facetVars <- facetVars[facetVars != "."]  # Sacar el puntito de las variables
            facetVars <- rev(facetVars)  # Dar vuelta, porque necesito primero los X y después los Y.

            # Conseguir los valores de las variables del facet usando información del brush.
            # brush_names <- c("panelvar1", "panelvar2", "nico", "panelvar3")
            brush_names <- names(input$scatterplot_brush)     # Get names
            panelvars_index <- grep("panelvar", brush_names)  # Get indexes of the panelvars
            panelvars <- brush_names[panelvars_index]         # Get the panelvars
            panelvals <- input$scatterplot_brush[panelvars]   # Get the values of the panelvars


            # Build the subset condition
            panelvals <- paste("'", panelvals, "'", sep = "") # Quote them, for later use in subset() with eval(parse())
            subset_condition = paste(paste(facetVars, "==", panelvals), collapse = " & ")

            # Subset the dataframe
            print(paste("-- Subsetting by facet condition:", subset_condition))
            d <- subset(d, eval(parse(text=subset_condition)))
          } else {
            print("-- Brush by facet mode OFF")
          }


          # Output an image if filtering returns a non-empty selection
          if(nrow(d) > 0) {
            print("-- Selection not empty: magick!")
            magick.cell <-  magickCell(d, p, ch=input$ch, 
                                       sortVar = input$x, 
                                       seed = values$seed, 
                                       n = n_max, 
                                       equalize_images = input$equalize_pics,
                                       normalize_images = input$normalize_pics,
                                       boxSize = boxSize, return_ucid_df = T)
            tmpimage <- magick.cell$img
            print(magick.cell$ucids)
          } else {
            # Output white if selection is empty
            print("-- Selection is empty")
            tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
          }

        } else {
          # Si no hay algo seleccionado
          print("-- Selection is empty")
          tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
        }

        # Esto es para que shiny pueda mostrar output de imagen
        # Ver si se puede pasar a TIFF
        tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
        list(src = tmpfile, contentType = "image/jpeg")
      }, deleteFile=TRUE)




      # Reactive image 2 ----------------
      # Reads input: input$position, input$x, input$y
      # Reads reactive: rv$pgnpts, values$cdata, values$seed
      # Writes reactive:
      output$pics2 <- renderImage({
        print("Rendering image 2: polygon selection")
        pgnpts <- rv$pgnpts

        # Si hay algo seleccionado
        if(nrow(pgnpts) >= 3){
          # Me quedo con las posiciones que me interesan
          # positions <- getPositions(input$position, max(paths$pos))
          positions <- rangeExpand(input$position, max(paths$pos))
          p <- subset(paths, pos %in% positions)
          d <- subset(values$cdata, pos %in% positions)
          if(!input$suspend_filters) d <- subset(d, filter == TRUE)

          # # Copiado del bloque del brush
          # if(input$facet != ""){
          #     # Si hay algo escrito en el facet formula field, filtrar el dataframe para mostrar solo las imagenes del facet seleccionado.
          #     facetVars <- all.vars(eval(parse(text=input$facet)))
          #     facetVars <- facetVars[facetVars != "."]  # Sacar el puntito de las variables
          #     facetVars <- rev(facetVars)  # Dar vuelta, porque necesito primero los X y después los Y.
          #
          #     # Conseguir los valores de las variables del facet usando información del brush.
          #     click_names <- names(input$vertex1)     # Get names
          #     panelvars_index <- grep("panelvar", click_names)  # Get indexes of the panelvars
          #     panelvars <- click_names[panelvars_index]         # Get the panelvars
          #     panelvals <- input$vertex1[panelvars]   # Get the values of the panelvars
          #
          #     print(names(input$vertex1))
          #     print(input$vertex1$panelvar1)
          #
          #     # Build the subset condition
          #     panelvals <- paste("'", panelvals, "'", sep = "") # Quote them, for later use in subset() with eval(parse())
          #     subset_condition = paste(paste(facetVars, "==", panelvals), collapse = " & ")
          #
          #     # Subset the dataframe
          #     d <- subset(d, eval(parse(text=subset_condition)))
          # }

          # Get relevant variables
          x <- input$x
          y <- input$y
          pts <- d[,c(x, y)]

          # Get points
          pips <- pip(pts, pgnpts)
          d <- d[as.logical(pips),]
          magick.cell <- magickCell(d, p, ch=input$ch, 
                                    sortVar = input$x, 
                                    seed = values$seed, 
                                    n = n_max, 
                                    equalize_images = input$equalize_pics,
                                    normalize_images = input$normalize_pics,
                                    boxSize = boxSize, return_ucid_df = T)
          tmpimage <- magick.cell$img
        } else {
          # Si no hay algo seleccionado output white
          tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
        }

        # Esto es para que shiny pueda mostrar output de imagen
        tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
        list(src = tmpfile, contentType = "image/jpeg")

        # TO-DO ####
        # use "magick.cell$ucids" to plot points for each sampled cells, overlaid in the main plot
      }, deleteFile=TRUE)



      # Reactive text 1 : hover ----------------
      # Reads input: input$hover$x, input$hover$y
      # Reads reactive:
      # Writes reactive:
      output$info <- renderText({
        tryCatch(
          expr = {
            paste0("x=", signif(input$hover$x, 3), "\t y=", signif(input$hover$y, 3))
            },
          error = function(cond) return("Waiting for mouse hover...")
        )
      })
      
      
      # Reactive image 3: pics_nearby ----------------
      # Reads input: values$points_nearby, input$position, input$x, input$y
      # Reads reactive: rv$pgnpts, values$cdata, values$seed
      # Writes reactive:
      output$pics_nearby <- renderImage({
        
        # points_nearby <- values$points_nearby
        
        # Si hay hover no NULL
        if(!is.null(input$hover[1])){
          # Me quedo con las posiciones que me interesan
          positions <- rangeExpand(input$position, max(paths$pos))
          p <- subset(paths, pos %in% positions)
          d <- subset(values$cdata, pos %in% positions)
          if(!input$suspend_filters) d <- subset(d, filter == TRUE)
          
          # El nearest hay que hacerlo despues de filtrar por facet
          print("-- Hover names:")
          print(names(input$hover))
          
          # Si hay algo escrito en el facet formula field, filtrar el dataframe para mostrar solo las imagenes del facet seleccionado.
          if(input$facet != "" && input$facet_brush){
            print("-- Brush by facet mode ON")
            
            # Cada "eje" del facet tiene los posibles valores de una variable de pdata.
            
            # Primero conseguir las variables en la fórmula (debería estar primera la variable en el "eje y del facet")
            facetVars <- all.vars(eval(parse(text=input$facet)))
            facetVars <- facetVars[facetVars != "."]  # Sacar el puntito de las variables
            facetVars <- rev(facetVars)  # Dar vuelta, porque necesito primero los X y después los Y.
            
            # Conseguir los valores de las variables del facet usando información del brush.
            # brush_names <- c("panelvar1", "panelvar2", "nico", "panelvar3")
            hover_names <- names(input$hover)     # Get names
            panelvars_index <- grep("panelvar", hover_names)  # Get indexes of the panelvars
            panelvars <- hover_names[panelvars_index]         # Get the panelvars
            panelvals <- input$hover[panelvars]   # Get the values of the panelvars
            
            
            # Build the subset condition
            panelvals <- paste("'", panelvals, "'", sep = "") # Quote them, for later use in subset() with eval(parse())
            subset_condition = paste(paste(facetVars, "==", panelvals), collapse = " & ")
            
            # Subset the dataframe
            print(paste("-- Subsetting by facet condition:", subset_condition))
            d <- subset(d, eval(parse(text=subset_condition)))
            
          } else {
            print("-- Brush by facet mode OFF")
          }
          
          
          # Find the closest point
          # d <- hover_closest(ui_input = input, cdata = d)
          print("-- Finding one of nearPoints, probably always by facet")
          d <- shiny::nearPoints(d, input$hover, maxpoints = 1, xvar = input$x, yvar = input$y)
          
          # Output an image if filtering returns a non-empty selection
          if(nrow(d) > 0) {
            print("-- Selection not empty: magick!")
            magick.cell <-  magickCell(d, p, ch=input$ch, 
                                       sortVar = input$x, 
                                       seed = values$seed, 
                                       n = 1, 
                                       equalize_images = input$equalize_pics,
                                       normalize_images = input$normalize_pics,
                                       boxSize = boxSize, return_ucid_df = T)
            tmpimage <- magick.cell$img
            print(magick.cell$ucids)
          } else {
            # Output white if selection is empty
            print("-- Selection is empty")
            tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
          }
          
        } else {
          # Si no hay algo seleccionado
          print("-- Selection is empty")
          tmpimage <- magick::image_blank(100,10,color = "white") %>% image_annotate(text = "Empty set")
        }
        
        # Esto es para que shiny pueda mostrar output de imagen
        # Ver si se puede pasar a TIFF
        tmpfile <- magick::image_write(tmpimage, tempfile(fileext='jpg'), format = 'jpg')
        list(src = tmpfile, contentType = "image/jpeg")
        
      }, deleteFile=TRUE)



      # Reactive table 1: FILTERS ----------------
      # Reads input: input$facet
      # Reads reactive: values$cdata
      # Writes reactive:
      output$filter_summary <- renderFormattable({

        text <- input$facet  # text = "~mpg+cyl"
        facetVars2 <- all.vars(eval(parse(text=text)))
        facetVars2 <- facetVars2[facetVars2 != "."]  # Sacar el puntito de las variables
        # facetVars2 <- paste(facetVars2, sep = ",")

        bind_rows(
          values$cdata %>%
            filter(filter) %>%
            select(!!facetVars2) %>%  # ok
            group_by(.dots = facetVars2) %>%
            summarise(count = n()) %>% mutate(stage = "filtered"),
          cdata %>%
            select(!!facetVars2) %>%  # ok
            group_by(.dots = facetVars2) %>%
            summarise(count = n()) %>% mutate(stage = "original")
        ) %>%
          tidyr::pivot_wider(id_cols = facetVars2,
                             names_from = stage,
                             values_from = count) %>%
          mutate(delta = round(filtered/original, digits = 2)) %>%
          # formattable(list(delta = formattable::color_tile("orange", "white")))
          formattable::formattable(list(
            delta = formatter("span", style = x ~ dplyr::case_when(
              x < .5 & x > .1 ~ style(display = "block", padding = "0 4px", `border-radius` = "4px", `background-color` = "orange"),
              x <= .1 ~ style(display = "block", padding = "0 4px", `border-radius` = "4px", `background-color` = "red"),
              TRUE ~ NA_character_)),
            original = formattable::color_tile("orange", "white")
          ))

        # http://enhancedatascience.com/2017/07/10/the-packages-you-need-for-your-r-shiny-application/
        # https://github.com/renkun-ken/formattable
      })
    }



    ### OBSERVADORES DE EVENTOS ###
    {
      # CLICK OBSERVER 1  ----------------
      observeEvent(
        # Observe clicks and add them to the polygon dataframe
        eventExpr = input$vertex1,  ## Double clicks
        handlerExpr = {
          writeLines("\nDoublelick event: Polygon point added")
          # Isolation may be futile in "observeEvent" contexts.
          
          # https://stackoverflow.com/questions/30588472/is-it-possible-to-clear-the-brushed-area-of-a-plot-in-shiny
          session$resetBrush("scatterplot_brush")  # https://shiny.rstudio.com/reference/shiny/0.14/session.html
          rv$brush_limits <- c(NULL, NULL, NULL, NULL)  # Trigger brush image update, so that it clears up.
          
          # Load all polygons
          pgnpts <- isolate(rv$pgnpts)  # will fire on click, isolation is prudent

          # Add a new row with the new point
          if(nrow(pgnpts) == 0){
            pgnpts <- data.frame(x = isolate(input$vertex1$x),
                                 y = isolate(input$vertex1$y),
                                 xvar = isolate(input$x),
                                 yvar = isolate(input$y),
                                 stringsAsFactors = F)
          } else {
            pgnpts <- bind_rows(pgnpts,
                                data.frame(x = isolate(input$vertex1$x),
                                           y = isolate(input$vertex1$y),
                                           xvar = isolate(input$x),
                                           yvar = isolate(input$y),
                                           stringsAsFactors = F))
          }

          rv$pgnpts <- pgnpts


        }, label = "Click observer 1")

      # CLICK OBSERVER 2  ----------------
      observeEvent(
        # Observe clicks and add them to the polygon dataframe
        eventExpr = input$vertex2,
        handlerExpr = {
          writeLines("\nClick event: Polygon cleared / Brush drawn")
          
          rv$pgnpts <- pgnpts_empty

        }, label = "Click observer 2")

      # BUTTON OBSERVER 1: FILTER ADD  ----------------
      observeEvent(
        # Acción al apretar el botón de agregar filtros
        eventExpr = input$add_filter,
        handlerExpr = {
          writeLines("\nAdd filter fired")
          # Isolation may be futile in "observeEvent" contexts.

          # Density 1D should be treated differently
          plot.type <- isolate(input$ptype)

          # Brush
          brush <- isolate(input$scatterplot_brush)
          if(!is.null(brush$xmin)){
            print("-- Brush mode ON")
            brpts <- square(brush$xmin, brush$ymin, brush$xmax, brush$ymax)
            brpts <- data.frame(x = brpts$x,
                                y = brpts$y,
                                xvar = isolate(input$x),
                                yvar = isolate(input$y),
                                type = input$filter_type,
                                stringsAsFactors = F)
            pgn <- brpts

            # Append the dataframe to the filters list
            filter_length <- length(rv$filters) + 1
            rv$filters[[filter_length]] <- pgn  # This will fire FILTER TAB EVENTS OBSERVER
          } else {
            print("-- Brush empty, no filter created from brush.")
          }

          # Polygon
          pgnpts <- rv$pgnpts
          if(nrow(pgnpts) >= 3){
            print("-- Polygon mode ON")
            pgn <- cbind(pgnpts, data.frame(type = input$filter_type, stringsAsFactors = F))

            # Append the dataframe to the filters list
            filter_length <- length(rv$filters) + 1
            rv$filters[[filter_length]] <- pgn  # This will fire FILTER TAB EVENTS OBSERVER
          } else {
            print("-- Not enough points for polygon filter.")
          }

          if (nrow(pgnpts) < 3 & is.null(brush$xmin)) {
            print("Incomplete or missing selection, ignoring button press.")
          } else {
            if(length(rv$filters) == 0){
              # This will fire FILTER TAB EVENTS OBSERVER
              updateCheckboxGroupInput(session = session, inputId = "stringFilters",
                                       choiceNames = c(),
                                       choiceValues = c()
              )
            } else {
              # This will fire FILTER TAB EVENTS OBSERVER
              updateCheckboxGroupInput(session = session, inputId = "stringFilters",
                                       choiceNames = paste0("polygon", 1:length(rv$filters)),
                                       choiceValues = 1:length(rv$filters), selected = c(as.numeric(input$stringFilters),filter_length)
              )
            }
          }
        })

      # BUTTON OBSERVER 2: EXIT  ----------------
      observeEvent(
      # Acción al apretar el botón de cerrar la app
        eventExpr = input[["quit"]],
        handlerExpr = {
          writeLines("\nQuit event fired!")

          saved <- list(cdata = values$cdata,
                        filters = rv$filters)

          # https://stackoverflow.com/questions/27365575/how-to-exit-a-shiny-app-and-return-a-value
          stopApp(saved)
        }
      )

    }


  }

