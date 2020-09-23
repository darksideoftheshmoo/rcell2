# Added to remove NOTES from devtools:check()
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# globalVariables(strsplit("channel choice counts h i level pos saved treatment w x y", " ")[[1]])
# "A horrible hack" (?)

# I have found the cause of re-rendering plots.
# clientData$output_scatterplot_width changes 10 PIXELS for no reason and triggers a re-render.
# Maybe its the scrollbar or something...

#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param cdata A Rcell data.frame (not the object).
#' @param pdata A "pdata" data.frame with position metadata.
#' @param paths Paths a la imagen de cada posición.
#' @param filters Vector de strings con los filtros. c() by default.
#' @param plotType "Hex", "Density", and "Dots" are available.
#' @param seed Seed value for sampling of cell images.
#' @param initial_facet Initial facet formula as a string.
#' @param facet_grid_option Use facet_grid (TRUE, default) or facet_wrap.
#' @param facets_scale_free Use facets with fixed scales (NULL, default) or free scales ("free").
#' @param ... Further arguments passed to magickCell()
#' @return Lots of stuff.
# @examples
# saved_data <- shinyCell(cdata, pdata, paths, plotType = "Dots")
#' @import shiny ggplot2 magick formattable shinydashboard
#' @importFrom grDevices rgb
#' @export
shinyCell <- function(cdata,
                      pdata,
                      paths,
                      filters = list(), filters.init_selected = T,
                      plotType = "Hex",
                      seed = 1,
                      initial_facet = "",
                      facet_grid_option = TRUE, facets_scale_free = NULL,
                      n_max = 100,
                      ...){

    # To-do
    # Invalid input$facet generates warnings and errors, this should be handled. Also, only "~", "." and "+" are handled in forumlas.
    # Implement more-than-2 variable faceting. The third and ith faceting variables of the brush are stored in "panelvar3" and so on (?)
    # Integrate polygon filter functionality, currently the drawn polygons do nothing (except show up).

    # runApp inicia la app inmediatamente, shinyApp solo no se dispara dentro de una función parece
    # saved <- runApp(shinyApp(ui, server))

    # server functions in a package .R script have an "enclosing environment" different from this function's (shinyCell) local environment.
    # that is one reason why shinyCell's argument dont reach the server function's env
    # another reason is that, somehow, shinyCell's env cannot be reached from within "runApp"
    # by replacing the server function's enclosing environment by shinyCell's env everything is fixed
    environment(shinyAppServer) <- environment() # https://stackoverflow.com/questions/44427752/distinct-enclosing-environment-function-environment-etc-in-r

    # shinyAppUI() is also defined as a function that returns a "fluidPage" object,
    # and thus suffers from the same problem as shinyAppServer()
    environment(shinyAppUI) <- environment()

    # Taken from example at Rbloggers
    # https://github.com/MangoTheCat/shinyAppDemo/blob/master/R/launchApp.R
    # Here shinyAppUI() must be executed in order to pass a fluidPage object to shinyApp/runApp
    saved <- shiny::runApp(list(ui = shinyAppUI(), server = shinyAppServer))

    # Imprimir cosas antes de cerrar la app
    print("Chau")

    # Devolver una lista con los objetos cdata cfilter y los stringFilters
    return(saved)
}
