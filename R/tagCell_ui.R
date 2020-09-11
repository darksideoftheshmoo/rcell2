#' Shiny app UI object for tagCell
#'
#' @import shiny shinyjs formattable shinydashboard
#'
tagCellUi <- function(){shiny::fluidPage(useShinyjs(),  # Set up shinyjs
  fluidPage(shiny::fluidRow(
    # shiny::column(
    #   width = 3,
    #   # shiny::p("A rough app to tag cells from CellID data"),
    #   shiny::selectInput('image_channel','Channel:', sort(unique(paths$channel)), "BF.out"),
    #   shiny::hr(),
    #   # shiny::verbatimTextOutput("info"),
    #   # shiny::verbatimTextOutput("logs"),
    #   shiny::p(shiny::actionButton(inputId = "quit", label = "Save & Quit"))
    # ),
    shiny::column(width = 9, 
                  offset = 0,
                  shiny::tabsetPanel(
                    # tags$head(tags$style("#{height:100vh !important;}")),
                    shiny::tabPanel("Cell pics",
                                    shiny::p("Tag the current cell to tag to the next one:"),
                                    shiny::p(shiny::plotOutput(outputId = "pics",
                                                               height = "100%", width = "100%"
                                                               # height = "auto", width = "auto"
                                    )),
                                    shiny::p(
                                      shiny::plotOutput(outputId = "plot"#, 
                                                        # height = "100%", width = "100%"
                                                        # height = "auto", width = "auto"
                                    )),
                                    shiny::p(
                                      shiny::tableOutput("saved_annotations")
                                    )
                                    ),
                    shiny::tabPanel("My ggplot")
                    )
                  ),
    shiny::column(width = 3, 
                  offset = 0,
                  shiny::p(
                    shiny::verbatimTextOutput("cell_ith"),
                    uiOutput("moreControls"),
                    shiny::actionButton(inputId = "prev_cell", label = "Previous", icon = shiny::icon("step-backward")),
                    shiny::actionButton(inputId = "next_cell", label = "Next", icon = shiny::icon("step-forward"))
                    ),
                  shiny::p(
                    shiny::hr(),
                    # shiny::selectInput('image_channel','Channel:', sort(unique(paths$channel)), "BF.out"),
                    shiny::actionButton(inputId = "prev_ucid", label = "Unskip", icon = shiny::icon("fast-backward")),
                    shiny::actionButton(inputId = "next_ucid", label = "Skip", icon = shiny::icon("fast-forward")),
                    shiny::hr(),
                    shiny::p(shiny::actionButton(inputId = "quit", label = "Save & Quit"))
                  ))
)))}
