#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param user_plot A ggplot object, ready to render.
#' @param debug_messages Print internal stuff, useful for debugging.
#' @param print_plot_on_exit Prints user_plot replacing it's data with the "brushed" subset
#' @return A data.frame with the brushed points from the original plot's data.
#' @import shiny tidyverse
#' @export
plotApp <- function(user_plot, debug_messages = F, print_plot_on_exit = F){

    environment(plotAppServer) <- environment()
    environment(plotAppUI) <- environment()
    saved <- shiny::runApp(list(ui = plotAppUI(), server = plotAppServer))
    
    if(print_plot_on_exit) print(user_plot %+% saved)

    # Devolver una lista con los objetos cdata cfilter y los stringFilters
    return(saved)
}
