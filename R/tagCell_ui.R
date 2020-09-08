#' Shiny app UI object for tagCell
#'
#' @import shiny formattable shinydashboard
#'
tagCellUi <- function(){shiny::fluidPage(
  fluidPage(
  shiny::fluidRow(
    shiny::column(width = 4,
           shiny::p("A rough app to tag cells from CellID data"),
           shiny::selectInput('image_channel','Channel:', sort(unique(paths$channel)), "BF.out"),
           shiny::hr(),

           shiny::verbatimTextOutput("info"),
           shiny::verbatimTextOutput("logs"),
           shiny::p(shiny::actionButton(inputId = "quit", label = "Save & Quit"))
    ),
    shiny::column(width = 8, offset = 0,
                  shiny::tabsetPanel(
                    # tags$head(tags$style("#{height:100vh !important;}")),
                    shiny::tabPanel("Cell pics",
                                    shiny::p("Tag the current cell to tag to the next one:"),
                                    shiny::p(shiny::plotOutput(outputId = "pics",
                                                               height = "100%", width = "100%"
                                                               # height = "auto", width = "auto"
                                    )),
                                    shiny::p(uiOutput("moreControls"),
                                             shiny::actionButton(inputId = "prev_cell", label = "Previous"),
                                             shiny::actionButton(inputId = "next_cell", label = "Next"),
                                             shiny::verbatimTextOutput("cell_ith"))
                    ),
                    shiny::tabPanel("My ggplot",
                                    shiny::plotOutput(outputId = "plot"#, 
                                                      # height = "100%", width = "100%"
                                                      # height = "auto", width = "auto"
                                    )
                    )
                  )
    )
  )
  )
)}
