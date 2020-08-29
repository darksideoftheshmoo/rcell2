# Added to remove NOTES from devtools:check()
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
globalVariables(strsplit("channel choice counts h i level pos saved treatment w x y", " ")[[1]])
# "A horrible hack" (?)

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
#' @import shiny ggplot2 magick
#' @importFrom grDevices rgb
#' @importFrom utils head
#' @export
shinyTag <- function(cdata,
                     pdata,
                     paths,
                     filters = list(),
                     plotType = "Hex",
                     seed = 1,
                     initial_facet = "",
                     facet_grid_option = TRUE, facets_scale_free = NULL,
                     ...){
  
  # To-do
  # Invalid input$facet generates warnings and errors, this should be handled. Also, only "~", "." and "+" are handled in forumlas.
  # Implement more-than-2 variable faceting. The third and ith faceting variables of the brush are stored in "panelvar3" and so on (?)
  # Integrate polygon filter functionality, currently the drawn polygons do nothing (except show up).
  
  # Load definitions
  source("R/definitions.R", local = T)
  
  #### UI ####
  source("R/tagCell_ui.R", local = T)
  #### SERVER ####
  source("R/tagCell_server.R", local = T)
  
  #### RUN APP ####
  # runApp inicia la app inmediatamente, shinyApp solo no se dispara dentro de una función parece
  # saved <- runApp(shinyApp(ui, server))
  shiny::runApp(list(ui=ui, server=server), launch.browser = TRUE)
  
  # Imprimir cosas antes de cerrar la app
  print("Chau")
  
  #### RETURN RESULT ####
  # Devolver una lista con los objetos cdata cfilter y los stringFilters
  return(saved)
}