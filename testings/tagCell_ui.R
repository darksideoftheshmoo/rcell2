#' Shiny app UI object for tagCell
#'
#' @import shiny formattable shinydashboard
#'
tagCellUi <- shiny::fluidPage(
  shiny::fluidRow(
    shiny::column(width = 4,
           shiny::p("A rough app to tag cells from CellID data"),
           shiny::selectInput('ch','Channel:', sort(unique(paths$channel)), "BF.out"),
           shiny::hr(),

           shiny::verbatimTextOutput("info"),
           shiny::verbatimTextOutput("logs"),
           shiny::p(shiny::actionButton(inputId = "quit", label = "Save & Quit"))
    ),
    shiny::column(width = 8, offset = 0,
                  shiny::tabsetPanel(
                    shiny::tabPanel("Cell pics",
                                    shiny::p("Tag the current cell to tag to the next one:"),
                                    shiny::plotOutput("pics", 
                                                      height = "100%", 
                                                      width = "100%")
                    )
                  )
    )
  )
)
