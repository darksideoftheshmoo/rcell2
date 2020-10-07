#' Shiny app UI object for shinyCell
#'
#' @import shiny formattable shinydashboard
#'
shinyAppUI <- function(){
  fluidPage(
    fluidRow(
      column(width = 4,
             p("A rough app to explore CellID data"),
             textInput(inputId = "position", label = "Positions:", value = paste(sort(unique(paths$pos)), collapse = ", ")),
             p(paste("Available positions are:",
                     numbers_to_intervals(pdata$pos),
                     " (leave blank for all)."
             )),
             # selectInput('ch','Channel:', sort(unique(paths$channel)), "BF.out"),
             selectInput('x','X variable', names(cdata)[!names(cdata) %in% names(pdata)],initial_vars[1]),  # Exclude pdata variables for plotting. Note: pos may be useful!
             selectInput('y','Y variable', names(cdata)[!names(cdata) %in% names(pdata)],initial_vars[2]),
             hr(),

             p(
              selectInput(inputId = "filter_type", choices = c("Additive", "Subtractive"),
                          label = "Filter type",
                          width = "100%"),
              actionButton(inputId = "add_filter", label = "Add filter"),
              actionButton(inputId = "quit", label = "Save and Quit"),
              tags$head(
                tags$style(HTML('#quit{background-color:orange}'))
              )
             ),

             verbatimTextOutput("info"),
             plotOutput("pics_nearby", height = "100%", width = "100%"),
             verbatimTextOutput("logs"),
             p()
      ),
      column(width = 8, offset = 0,
             plotOutput(outputId = "scatterplot",
                        brush = brushOpts(id = "scatterplot_brush", fill = "#ccc", resetOnNew = TRUE),
                        click = "vertex2",
                        dblclick = "vertex1",
                        hover = "hover",
                        height = "100%",
                        width = "100%"
                        #inline = TRUE
             )
      )
    ),

    hr(),

    column(width = 12,
           tabsetPanel(
             # tabPanel("Brush pics",
             tabPanel("Brush and Polygon pics",
                      column(width = 2,
                             p("Image options:"),
                             selectInput('ch','Channel:', sort(unique(paths$channel)), "BF.out"),
                             checkboxInput(inputId = "normalize_pics", label = "Normalize", value = FALSE),
                             checkboxInput(inputId = "equalize_pics", label = "Equalize", value = FALSE)
                      ),
                      column(width = 10,
                             p("Cells in brushed points:"),
                             plotOutput("pics", height = "100%", width = "100%"),  # Jugando un poco con los tamanos, todavia no sale lindo
                             plotOutput("pics2", height = "100%", width = "100%")  # Jugando un poco con los tamanos, todavia no sale lindo
                      )
             ),
             # tabPanel("Poly pics",
             #          p("Cells in polygon:"),
             #          plotOutput("pics2", height = "100%", width = "100%")  # Jugando un poco con los tamanos, todavia no sale lindo
             # ),
             
             # tabPanel("Image settings",
             #          column(width = 4,
             #                 checkboxInput(inputId = "normalize_pics", label = "Normalize images", value = FALSE),
             #                 checkboxInput(inputId = "equalize_pics", label = "Equalize images", value = FALSE)
             #                 # shiny::numericInput("n_images_max", "Max. number of cells per tile", 10)
             #          )
             # ),

             tabPanel("Filter Settings",
                      column(width = 4,
                      checkboxInput(inputId = "overlay_polygons", label = "Draw filters?", value = TRUE),

                      checkboxInput(inputId = "suspend_filters", label = "Suspend filters?", value = FALSE),

                      selectInput("truth_mode", "Filter mode", c("Discard > Keep" = "all", "Keep > Discard" = "any"), selected = "all"),
                      # selectInput("cell_unique_id_field", "Filter primary key", c("ucid" = "ucid", "ucid_time" = "ucid_time"), selected = "all"),


                      #uiOutput("filters")
                      # Fue necesario porque la lazy evaluation no cargaba los filtros hasta que abriera este tab por primera vez

                      # Using "names(filters)" as temporary fix, while switching from argument class c() -> list()
                      checkboxGroupInput("stringFilters", "Choose filters to apply:",
                                         # choiceNames= length(filters),
                                         # choiceValues= length(filters),
                                         selected = if(filters.init_selected) 1:length(filters) else c(),
                                         choiceValues = if(length(filters)==0) c() else 1:length(filters),
                                         choiceNames = if(length(filters)==0) c() else paste("polygon", 1:length(filters)),
                                         width = '100%')
                      ),
                      column(width = 8,
                             # tableOutput('filter_summary')
                             formattableOutput("filter_summary")
                             )
             ),

             tabPanel("Plot Settings",
                      # selectInput("facets_scale_free", "Facet scale", list(fixed=NULL, free="free"), NULL),
                      column(width = 6,
                             sliderInput("plotDimX",label = "Plot X dimentions",
                                         min = 500, max = 2000, step = 50, width = '100%', post = " px", value = 600),
                             sliderInput("plotDimY",label = "Plot Y dimentions",
                                         min = 500, max = 2000, step = 50, width = '100%', post = " px", value = 600)),
                      column(width = 6,
                             selectInput("ptype", "Plot type", c("Hex", "Density", "Dots"), plotType),
                             textInput(inputId = "facet", label = "Facet formula:", value = initial_facet),
                             p(paste("Use up to two variables:", paste(names(pdata), collapse = ", ")), "."),  # "pos" may be nice for faceting
                             p("Use a dot for one-variable facets. Leave blank for none."),
                             checkboxInput(inputId = "facet_brush", label = "Brush by facet?", value = TRUE)
                      )
             )
           )
    )
  )
}
