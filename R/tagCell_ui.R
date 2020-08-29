#' Shiny app UI object for tagCell
#'
#' @import shiny formattable shinydashboard
#'
tagCellUi <- shiny::fluidPage(
  shiny::fluidRow(
    shiny::column(width = 4,
           shiny::p("A rough app to tagg cells from CellID data"),
           shiny::textInput(inputId = "facet", label = "Facet formula:", value = ""),
           shiny::checkboxInput(inputId = "facet_brush", label = "¿Brush by facet?", value = TRUE),
           shiny::p(paste("Use a dot for one-variable facets. Leave blank for none.",
                   "Use up to two variables:",
                   paste(names(pdata), collapse = ", ")), "."),  # "pos" may be nice for faceting
           shiny::textInput(inputId = "position", label = "Positions:", value = paste(sort(unique(paths$pos)), collapse = ", ")),
           shiny::p(paste("Available positions are:",
                   paste(sort(pdata$pos), collapse = ", "),
                   ". Leave blank for all."
           )),
           shiny::selectInput('ch','Channel:', sort(unique(paths$channel)), "BF.out"),
           shiny::selectInput('x','X variable', names(cdata)[!names(cdata) %in% names(pdata)],"a.tot"),  # Exclude pdata variables for plotting. Note: pos may be useful!
           shiny::selectInput('y','Y variable', names(cdata)[!names(cdata) %in% names(pdata)],"fft.stat"),
           shiny::hr(),
           shiny::p("Filter actions:", shiny::actionButton(inputId = "add_filter", label = "Add filter")),
           shiny::p(shiny::selectInput(inputId = "filter_type", choices = c("Additive", "Subtractive"), label = "Filter type")),

           shiny::verbatimTextOutput("info"),
           shiny::verbatimTextOutput("logs"),
           shiny::p(shiny::actionButton(inputId = "quit", label = "Save & Quit"))
    ),
    shiny::column(width = 8, offset = 0,
                  shiny::tabsetPanel(
                    shiny::tabPanel("Brush pics",
                                    shiny::p("Cells in brushed points:"),
                                    shiny::plotOutput("pics", 
                                                      height = "100%", 
                                                      width = "100%")  # Jugando un poco con los tamaños, todavía no sale lindo
                    )
                  )
    )
  )
)
