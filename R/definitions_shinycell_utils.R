#' Plot shinyCell polygon filters
#' 
#' Useful to check out what areas the filters are covering.
#'
#' @param saved_data The output of shinyCell.
# @param .type a string, either "Subtractive" or "Additive", the two types of filters.
#' @param print_plots Set to false to prevent printing the plots on execution.
#'
#' @importFrom rlang parse_expr
#'
#' @return Prints a lis
#' @export
#'
# @examples
plot_filters <- function(saved_data,
                         # .type = "Subtractive",
                         print_plots = TRUE){
  
  pgnfilters.o <- saved_data$filters %>% bind_rows(.id = "polygon") # %>% filter(type == .type)
  
  pgn.vars <- pgnfilters.o %>% select(xvar, yvar) %>% unique() %>% bind_rows(data.frame(xvar = "el.p", yvar="a.tot"))
  variables <- pgn.vars %>% select(xvar, yvar) %>% plyr::adply(.margins = 1, function(x){
    asd <- c(xvar=x$xvar, yvar=x$yvar)
    asd[order(asd)]
  }) %>% unique()
  
  plot_list <- list()
  
  for(i in 1:nrow(variables))local({
    .x <- variables[i,1]
    .y <- variables[i,2]
    x_ <- rlang::parse_expr(.x)
    y_ <- rlang::parse_expr(.y)
    
    
    pgnfilters.o.unswapped <- pgnfilters.o %>% 
      filter(xvar == .x & yvar == .y) %>% 
      dplyr::rename(!!x_ := x,
                    !!y_ := y)
    pgnfilters.o.swapped <- pgnfilters.o %>% 
      filter(xvar == .y & yvar == .x) %>% 
      dplyr::rename(!!x_ := y,
                    !!y_ := x)
    
    p <- bind_rows(pgnfilters.o.unswapped, pgnfilters.o.swapped) %>%
      ggplot() +
      geom_point(aes(
        x = !!x_,
        y = !!y_
      ), data = filter(saved_data$cdata, t.frame == 0)) +
      geom_polygon(aes(x = !!x_, 
                       y = !!y_,
                       color = polygon,
                       linetype = type),  # No se puede usar "fill", conflictua con "scale_fill" en Hex y Density
                   size = 1,
                   alpha = .1) +
      theme_minimal()
    
    if(print_plots) print(p)
    
    plot_list[[paste(.x,.y,sep="_")]] <<- p
  })
  
  return(plot_list)
}