# Added to remove NOTES from devtools:check()
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# globalVariables(strsplit("channel choice counts h i level pos saved treatment w x y", " ")[[1]])
# "A horrible hack" (?)

# I have found the cause of re-rendering plots.
# clientData$output_scatterplot_width changes 10 PIXELS for no reason and triggers a re-render.
# Maybe its the scrollbar or something...

#' Safe select column
#' 
#' Checks whether the column selection with \code{[[]]} is null before returning.
#' 
safe_select <- function(.df, .name){
  
  if(!all(.name %in% names(.df))){
    .name_missing <- .name[!.name %in% names(.df)]
    stop(paste("Error, selected columns'", .name_missing,  "'does not exist.", collapse = " "))
  }
  
  column <- .df[[.name]]
  
  if(is.null(column)) stop(paste0("Error, selected column '", .name,  "' does not exist."))
  
  return(column)
}

#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param cdata A Rcell "cdata" data.frame with the CellID variables (not the cell.data object, \code{cell.data$data}).
#' @param pdata An optional "pdata" data.frame, with positions' metadata (NULL by default).
#' @param paths A "paths" data.frame, with paths to each positions' images (i.e. \code{cell.data$images}).
#' @param filters An optional list with the filters from a previous shinyCell run (dataframes with points of 2D polygons). An empty \code{list()} by default.
#' @param plotType Type of the filtering plot, either: "Dots" (point scatterplot, defaut), "Hex" (2D histogram with hexagonal bins), "Density", "Pics" (a \code{cellSpread} plot).
#' @param seed Seed value for sampling of cell images.
#' @param initial_facet Initial ggplot facet formula as a string (for example: "~pos+group_1")
#' @param initial_vars Initial cdata variables as a character vector (defaults to \code{c('a.tot', 'fft.stat')}).
#' @param facet_grid_option Use ggplot's facet_grid (TRUE, default) or facet_wrap (FALSE).
#' @param facets_scale_free Use ggplot's facets with fixed scales (NULL, default) or free scales ("free").
#' @param boxSize Size in pixels for individual cells' images.
#' @param filter_progress_file Path to an RDS file, used for saving filtering progress (in case something goes wrong). Using FALSE disables this feature. Set to NULL (the default) to let tempfile() choose a path for the RDS, or set to a valid path of your choice.
#' @param launch.browser Set to \code{'firefox'} or equivalent to launch the app in-browser (\code{FALSE} by default). Useful when launching fails with error \code{Error in utils::browseURL(appUrl)} or similar.
#' @param skip_input_check If FALSE (default)
#' @param ... Further arguments passed to \code{magickCell()}.
#' @return A named list with the original cdata and a list of filters. The cdata includes an extra "filter" column, indicating if a row is to be kept (TRUE) or filtered out (FALSE). The list of filters can be passed as a filter argument, and can be plotted with \code{plot_filters}.
#' @examples
#' 
#' # Minimal example:
#' 
#' path <- "/path/to_your/cellid_images/"
#' 
#' cell.data <- rcell2::cell.load.alt(path = path)
#' 
#' cdata <- cell.data$data  # CellID dataframe
#' 
#' images <- cell.data$images  # Image paths
#' 
#' pdata <- read.csv("data/pdata.csv")  # "Position" metadata
#' 
#' filter.output <- 
#'   rcell2::shinyCell(cdata = cdata, 
#'                     pdata = pdata, 
#'                     paths = images)
#'                     
#' plot_filter(filter.output)
#' 
#' cdata.filtered <- dplyr::filter(filter.output, filter)
#'   
#' @import shiny ggplot2 magick formattable shinydashboard
#' @importFrom grDevices rgb
#' @export
shinyCell <- function(cdata,
                      pdata=NULL,
                      paths,
                      filters = list(),
                      filters.init_selected = T,
                      plotType = "Dots",
                      seed = 1,
                      initial_facet = "", 
                      initial_vars = c("a.tot", "fft.stat"),
                      facet_grid_option = TRUE,
                      facets_scale_free = "fixed",
                      n_max = 100, boxSize = 80,
                      filter_progress_file = NULL,
                      launch.browser = F,
                      skip_input_check = F,
                      ...){
  
  if(!skip_input_check){
    if(has.na_nan_inf(cdata)) stop("Error: your 'cdata' dataframe has NaN, NA and/or Inf values. Try using 'rcell2:::has.na_nan_inf() to find problematic columns.'")
    if(has.na_nan_inf(pdata)) stop("Error: your 'pdata' dataframe has NaN, NA and/or Inf values. Try using 'rcell2:::has.na_nan_inf() to find problematic columns.'")
    if(has.na_nan_inf(paths)) stop("Error: your 'paths' dataframe has NaN, NA and/or Inf values. Try using 'rcell2:::has.na_nan_inf() to find problematic columns.'")
  }
  
  if(is.null(pdata)) pdata <- data.frame(pos = safe_select(cdata, "pos"))
  
  if(!all(names(pdata) %in% names(cdata))) stop("Error: cdata does not contain names in pdata, join them first :)")
  if(!is.character(facets_scale_free)) stop("Error: facets_scale_free must be a string accepted by ggplot's scales argument. See ?facet_wrap.")
  if(!is.null(initial_vars)) {
      if(!all(initial_vars %in% names(cdata))) stop("Error: cdata does not contain some of the initial_vars")
      if(!is.character(initial_vars)) stop("Error: initial_vars is not a character vector")
      if(length(initial_vars) != 2) stop("Error: initial_vars must be of length 2 (for the horizontal and vertical axes).")
  }
  
  if(is.null(filter_progress_file)) {
    filter_progress_file <- tempfile(pattern = "shinyCell_progress", fileext = ".RDS")
    print(paste("-- Saving filter progress to temporary file:", filter_progress_file))
  }
    
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
  if(!isFALSE(launch.browser)) {
    default_browser <- getOption("browser")
    options(browser = launch.browser)
    
    saved <- shiny::runApp(list(ui = shinyAppUI(), 
                                server = shinyAppServer), 
                           launch.browser = T)
    options(browser = default_browser)
    
  } else {
    saved <- shiny::runApp(list(ui = shinyAppUI(), 
                                server = shinyAppServer))
  }


  # Imprimir cosas antes de cerrar la app
  print("Chau! returning 'invisible' results...")
  
  # Append progress file path
  if(!isFALSE(filter_progress_file)) saved$filter_progress_file <- filter_progress_file

  # Devolver una lista con los objetos cdata cfilter y los stringFilters
  return(invisible(saved))
}

