# Added to remove NOTES from devtools:check()
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# globalVariables(strsplit("channel choice counts h i level pos saved treatment w x y", " ")[[1]])
# "A horrible hack" (?)

# I have found the cause of re-rendering plots.
# clientData$output_scatterplot_width changes 10 PIXELS for no reason and triggers a re-render.
# Maybe its the scrollbar or something...

#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param cdata dataframe of "cell data"
#' @param pdata dataframe "position data"
#' @param paths dataframe of image paths 
#' @param cell_tags list of named vectors corresponding to tag groups and tags. 
#' @param tag_box_size size of the image crop in pixels
#' @param cell_resize resize of the image crop in pixels
#' @param tag_channels_select a vector giving names for the image channels: c("BF", "YFP.out", etc....)
#' @param n_max max number of boxes in the image
#' @param seed seed for random sampling of images
#' @param tmp_csv_output file path into which tagging information will be dumped progressively
#' @param tag_ggplot a ggplot object to display in the second tab, may be used for something someday.
# @param ... extra arguments, not used.
#' @return Lots of stuff.
# @examples
# saved_data <- shinyCell(cdata, pdata, paths, plotType = "Dots")
#' @import shiny ggplot2 magick
#' @importFrom grDevices rgb
#' @importFrom utils head
#' @export
tagCell <- function(cdata,
                    pdata,
                    paths,
                    cell_tags,
                    tag_box_size = 50,
                    cell_resize=100,
                    tag_channels_select=c("BF"),
                    n_max=10,
                    seed = 1,
                    tmp_csv_output=tempfile(tmpdir = "./", fileext = ".txt"),
                    tag_ggplot = NULL,
                    # plotType = "Hex",
                    # initial_facet = "",
                    # facet_grid_option = TRUE,
                    # facets_scale_free = NULL,
                    ...){
  
  # To-do
  # Invalid input$facet generates warnings and errors, this should be handled. Also, only "~", "." and "+" are handled in forumlas.
  # Implement more-than-2 variable faceting. The third and ith faceting variables of the brush are stored in "panelvar3" and so on (?)
  # Integrate polygon filter functionality, currently the drawn polygons do nothing (except show up).
  
  environment(tagCellServer) <- environment()
  environment(tagCellUi) <- environment()
  
  #### RUN APP ####
  # runApp inicia la app inmediatamente, shinyApp solo no se dispara dentro de una función parece
  # saved <- runApp(shinyApp(ui, server))
  saved <- shiny::runApp(list(ui = tagCellUi(), server = tagCellServer))
  
  # Imprimir cosas antes de cerrar la app
  print("Chau")
  
  #### RETURN RESULT ####
  # Devolver una lista con los objetos cdata cfilter y los stringFilters
  return(saved)
}