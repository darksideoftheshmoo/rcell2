#' Check dataframes for NA, NaN or Inf values.
#' 
#' @param df The data.frame to check.
#' @param print.which Print which columns have bad values to the console.
#' @return Logical value: TRUE (chekc passed) or FALSE (check failed; nasty values found).
#' @keywords internal
has.na_nan_inf <- function (df, print.which = F) {
  r <- lapply(df, function(x) is.nan(x) | is.infinite(x) | 
                is.na(x))
  result <- any(unlist(r))
  
  if(print.which & result){
    r.cols <- unlist(lapply(r, any))
    
    if(result){
      print("Found either NA, NaN or Inf values in the following columns:")
    } else {
      print("Did not find either NA, NaN or Inf values in the input's colums.")
    }
    
    print(r.cols[r.cols])
  }
  
  return(result)
}

#' Plot shinyCell polygon filters
#' 
#' Useful to check out what areas the filters are covering.
#'
#' @param saved_data The output of shinyCell, or a list: \code{list(cdata = NULL, filter = saved_data$filters)}. Note that the only effect of \code{cdata = NULL} is that points will not be drawn.
# @param .type a string, either "Subtractive" or "Additive", the two types of filters.
#' @param print_plots Set to false to prevent printing the plots on execution.
#'
#' @importFrom rlang parse_expr
#' @importFrom plyr adply
#'
#' @return A list of ggplots ready to print.
#' @export
#'
# @examples
plot_filters <- function(saved_data,
                         # .type = "Subtractive",
                         print_plots = TRUE){
  
  pgnfilters.o <- saved_data$filters %>% 
    bind_rows(.id = "polygon") # %>% filter(type == .type)
  
  pgn.vars <- pgnfilters.o %>% select(xvar, yvar) %>% unique() #%>% bind_rows(data.frame(xvar = "el.p", yvar="a.tot"))
  variables <- pgn.vars %>% select(xvar, yvar) %>% plyr::adply(.margins = 1, function(x){
    var_names <- c(xvar=x$xvar, yvar=x$yvar)
    var_names[order(var_names)]
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
    
    d <- bind_rows(pgnfilters.o.unswapped, pgnfilters.o.swapped)
    
    p <- d %>%
      ggplot() +
      # geom_point(aes(x = !!x_,y = !!y_), data = filter(saved_data$cdata, t.frame == 0)) +
      {if(!is.null(saved_data$cdata)) geom_point(aes(x = !!x_,y = !!y_), data = saved_data$cdata) else NULL} +
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

#' Bind shinyCell polygon filters by variable pairs
#' 
#' Unites polygons with matching dimensions. Useful to check out what areas the filters are covering (see \code{plot_bound_filters}).
#'
#' @param saved_data The output of shinyCell.
# @param .type a string, either "Subtractive" or "Additive", the two types of filters.
#'
#' @importFrom rlang parse_expr
#' @importFrom plyr adply
#'
#' @return A list of polygons bound by variable, with names unique to variable pairs (by sort). Note that "x" and "y" column names may be swapped relative to the input order.
#' @export
#'
# @examples
bind_filters <- function(saved_data,
                         # .type = "Subtractive",
                         print_plots = TRUE){
  
  pgnfilters.o <- saved_data$filters %>% 
    bind_rows(.id = "polygon") # %>% filter(type == .type)
  
  pgn.vars <- pgnfilters.o %>% select(xvar, yvar) %>% unique() #%>% bind_rows(data.frame(xvar = "el.p", yvar="a.tot"))
  variables <- pgn.vars %>% select(xvar, yvar) %>% plyr::adply(.margins = 1, function(x){
    var_names <- c(xvar=x$xvar, yvar=x$yvar)
    var_names[order(var_names)]
  }) %>% unique()
  
  filter_list <- list()
  
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
    
    d <- bind_rows(pgnfilters.o.unswapped, pgnfilters.o.swapped) %>% 
      dplyr::rename(x = !!x_,
                    y = !!y_)
    
    filter_list[[paste(.x,.y,sep="_")]] <<- d
  })
  
  return(filter_list)
}

#' Plot bound shinyCell polygon filters by variable pairs
#' 
#' Useful to check out what areas the filters are covering.
#'
#' @param bound_filters The output of \code{bind_filters}.
#'
#' @import ggplot2
#'
#' @return Plots for the polygons in bound_filters.
#' @export
#'
# @examples
plot_bound_filters <- function(bound_filters){
  plots <- 
    lapply(bound_filters, function(bfs){
      ggplot(bfs) +
        geom_polygon(aes(x=x,y=y,color=polygon,linetype=type), alpha = 0) + ggtitle(paste(bfs$xvar[1], bfs$yvar[1])) +
        theme_minimal()}
    )
  
  plots
}
