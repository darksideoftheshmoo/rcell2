#' Shiny app UI object for shinyCell
#'
#' @import shiny 
# @import DT
#'
plotAppUI <- function(){
  fluidPage(
    shiny::fluidRow(
      shiny::actionButton(inputId = "quit", label = "Save & Quit")
    ),
    shiny::fluidRow(#width = 6,
      plotOutput(outputId = "scatterplot",
                 brush = brushOpts(id = "scatterplot_brush", fill = "#ccc", resetOnNew = F),
                 click = "vertex2",
                 dblclick = "vertex1",
                 hover = "hover"
      )),
    shiny::fluidRow(#width = 6,
      shiny::dataTableOutput("selection_table")
    )
  )}
