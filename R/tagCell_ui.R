#' Shiny app UI object for tagCell
#'
#' @import shiny shinyjs formattable shinydashboard
#'
tagCellUi <- function(){shiny::fluidPage(shinyjs::useShinyjs(),  # Set up shinyjs
                                         tags$head(
                                           tags$style(HTML('#quit{background-color:orange}'))
                                         ),
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
                                    shiny::p("Tag the current cell using the left panels."),
                                    shiny::p(
                                      # Cell magick images
                                      shiny::plotOutput(outputId = "pics",
                                                        height = "100%",
                                                        width = "100%"
                                                        # height = "auto", width = "auto")
                                    )),
                                    shiny::p(verbatimTextOutput("hover_info")),
                                    shiny::p(
                                      # User plot
                                      shiny::plotOutput(outputId = "plot", 
                                                        # height = "100%", width = "100%"
                                                        # height = "auto", width = "auto"
                                                        height = "200px",
                                                        click = "plot_click",
                                                        hover = "plot_hover"),
                                      # Cell strips
                                      shiny::plotOutput(outputId = "pics2", 
                                      # shiny::plotOutput(outputId = "pics2",
                                                        # height = "100%",
                                                        # width = "100%"
                                      )
                                    )
                                    ),
                    shiny::tabPanel("Saved annotations",
                                    shiny::tableOutput("saved_annotations")
                                    )
                    )
                  ),
    shiny::column(width = 3, 
                  offset = 0,
                  shiny::p(
                    "Current cell and frame info:",
                    shiny::verbatimTextOutput("cell_ith"),
                    uiOutput("moreControls")
                    ),
                  shiny::p(shiny::hr()),
                  shiny::p("Change frames using these buttons:"),
                  shiny::p(
                    shiny::actionButton(inputId = "prev_cell", label = "Previous", icon = shiny::icon("step-backward")),
                    shiny::actionButton(inputId = "next_cell", label = "Next", icon = shiny::icon("step-forward")),
                  ),
                  shiny::p("Or click in the plot to jump between frames."),
                  shiny::p(shiny::hr()),
                  shiny::p("Jump to another cell using these:"),
                  shiny::p(
                    # shiny::selectInput('image_channel','Channel:', sort(unique(paths$channel)), "BF.out"),
                    shiny::actionButton(inputId = "prev_ucid", label = "Unskip", icon = shiny::icon("fast-backward")),
                    shiny::actionButton(inputId = "next_ucid", label = "Skip", icon = shiny::icon("fast-forward")),
                    shiny::hr(),
                    shiny::p(shiny::actionButton(inputId = "save", label = "Write progress")),
                    shiny::hr(),
                    shiny::p(shiny::actionButton(inputId = "quit", label = "Save & Quit"))
                  )
    )
)))}
