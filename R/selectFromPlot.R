#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param user_plot A ggplot
#' @param debug_messages print stuff
#' @return Lots of stuff.
#' @import shiny tidyverse
#' @export
plotApp <- function(user_plot, debug_messages = F){

    environment(plotAppServer) <- environment()
    environment(plotAppUI) <- environment()
    saved <- shiny::runApp(list(ui = plotAppUI(), server = plotAppServer))
    
    print(user_plot %+% saved)

    # Devolver una lista con los objetos cdata cfilter y los stringFilters
    return(saved)
}
