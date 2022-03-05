### DEFINITIONS ###

#' Find cell closest to hover point in shiny
#' @param ui_input the "input" object from whiny server
#' @param cdata the dataframe to filter
#' @keywords internal
hover_closest <- function(ui_input, cdata){
    d.closest <- sqrt(  
        (ui_input$hover$x - cdata[[ui_input$x]])^2 +
            (ui_input$hover$y - cdata[[ui_input$y]])^2
    ) %>% {cdata[which.min(.),]}
    return(d.closest)
}

#' convert sequence of numbers to ranges
#' @importFrom IRanges IRanges
#' @importFrom IRanges reduce
#' @keywords internal
numbers_to_intervals <- function(numbers = c(1:10, 20:30)){
    print("F0.1 numbers_to_intervals")
    x <- IRanges::reduce(IRanges::IRanges(numbers, width = 1))
    x <- as.data.frame(x)[,-3]
    x <- apply(x, 1, function(x) paste(x, collapse = "-"))
    paste(x, collapse = ", ")
}

#' rangeExpand
#' 
#' Function for producing an integer vector, by parsing integer ranges provided as a string: a comma-separated list of integers and "ranges" (i.e. "1:3" or "1-3").
#' 
#' See https://www.rosettacode.org/wiki/Range_expansion#R
#' @param text the input string to parse. Expect unexpected behaviour with decimal numbers, use integers only: \code{text = "1:2, 7-9"}.
#' @param maxPos total amount of positions (an int)
#' @return Integer vector, an expansion of the ranges specified in the input string.
#' @keywords internal
# @export
rangeExpand <- function(text = "1:2, 7-9", maxPos) {
    print("F0.2 rangeExpand")
    text <- gsub("[^0-9,\\:\\-]", replacement = "", text)  # Remove anything that is not integers or separators: ":" "-" ","
    lst <- gsub("(\\d)[-\\:]", "\\1:", unlist(strsplit(text, ",")))
    numeritos <- unlist(sapply(lst, function (x) eval(parse(text=x))), use.names=FALSE)
    # numeritos <- numeritos[numeritos < maxPos]  # better to fail
    if (length(numeritos) > 0) as.numeric(numeritos) else 1:maxPos
}

#' Return Shiny brush vertices
#' @keywords internal
square <- function(x1, y1, x2, y2){
    print("F0.3 square")
    return(list(
        x = c(x1, x1, x2, x2),
        y = c(y1, y2, y2, y1)
    ))
}

# To-do: hacer algo para poder armar filtros fuera de la app. ¿Cómo hacer filtros de una sola variable? El polígono no sirve.
#filterBox <- function(xvar, xmin, xmax, yvar, ymin, ymax)

#' Parse ggplot facet formula and get variables as a character vector
#' 
#' Una función para procesar el facet string y que devuelva las variables presentes en names(pdata) en un character vector
#' 
#' @keywords internal
getFacetVars <- function(pdata, facetFormulaString = "pos ~ treatment"){
    print("F1 getFacetVars")
    facetFormula <- eval(parse(text=facetFormulaString))
    facetVars <- all.vars(facetFormula)
    facetVars <- facetVars[facetVars %in% names(pdata)]
    return(facetVars)
}

#' Get a list of numbers from the "position" input string
#' 
#' @keywords internal
#' 
getPositions <- function(input, numPos){
    print("F2 getPositions")
    positions <- sort(strtoi(strsplit(input, split = ' |,|,,| ,|, ')[[1]]))
    if(length(positions) == 0) positions <- seq.int(from = 1, to = numPos, by = 1)
    return(positions)
}

#' Apply polygonal filters to cdata
#' 
#' Adds a boolean "filter" column to the "cdata" dataframe, based on a polygon list (typically output by shinyCell).
#' 
#' @param polygon_df_list A list of polygon dataframes with columns: x (values) y (values) xvar (variable name for x values) yvar (variable name for y values) type ("Subtractive" or "Additive")
#' @param cdata A "cdata" dataframe.
#' @param truthMode Logical priority for "Subtractive" and "Additive" polygon filter types, passed to \code{calculateTruth}. Should be one of "all" (Subtractive overcomes Additive) or "any" (Additive overcomes Subtractive).
#' @param cell_unique_id_field Name for the column holding the unique identifier (a "primary key") for each data point (i.e. the "ucid" is not suficcient for time series datasets).
#' 
#' @return a "saved_data" list object, where the cdata is appended a "filter" logical column.
#' 
polyFilterApply <- function(polygon_df_list,
                            cdata,
                            truthMode = "all",
                            cell_unique_id_field = "ucid"){
  
  # Check uniqueness of cell_unique_id_field
  is_primary_key <- length(cdata[,cell_unique_id_field, drop = T]) == length(unique(cdata[,cell_unique_id_field, drop = T]))
  if(!is_primary_key) stop(paste0("Error in polyFilterApply: '", cell_unique_id_field, "' is not a unique identifier."))

  print("F3 polyFilterApply")
  # Initialize empty cfilter, only with a primary key and TRUE filter
  cfilter <- data.frame(id = cdata[,cell_unique_id_field], 
                        filter = TRUE)
  names(cfilter)[1] <- cell_unique_id_field

  # Populate cfilter columns, append one column per filtering polygon
  print("F4.0 polyFilterApply")
  for (i in seq_along(polygon_df_list)) cfilter = do.call(what = polyFilterCell,
                                                          args = list(cdataDF = cdata,
                                                                      filterDF = cfilter,  # Iterates over the output
                                                                      polygonDF = polygon_df_list[[i]],
                                                                      polygonName = paste0("polygon", i)))
  
  # Recalcular la verdad
  print("F6.0 polyFilterApply")
  cfilter <- calculateTruth(filterDF = cfilter, mode = truthMode, cell_unique_id_field = cell_unique_id_field)

  # Add TRUTH column to cdata
  print("F7.0 polyFilterApply")
  cdata <- applyFilter(cdata, cfilter, 
                       cell_unique_id_field = cell_unique_id_field)

  # Return cdata and cfilter in a list.
  return(list(cdata = cdata,
              cfilter = cfilter))
}

#' Build filterDF from polygonDF and cdataDF
#' @keywords internal
polyFilterCell <- function(cdataDF, filterDF, polygonDF, polygonName = 1){
    print(paste("F4: polygon name:", polygonName))
    # polygonDF is NULL

    xvar <- polygonDF[1, "xvar", drop = T]
    yvar <- polygonDF[1, "yvar", drop = T]

    pips <- pip(points = cdataDF[,c(xvar, yvar)], pgn = polygonDF)

    type <- polygonDF[1, "type"]  # Poolygon filter type, as chosen at the shiny app

    if(type == "Additive") {
        filterDF[, paste0(polygonName, "_", type)] <-  as.logical(pips)

    } else if(type == "Subtractive") {
        filterDF[, paste0(polygonName, "_", type)] <- !as.logical(pips)

    } else {
        print("F4: polygon name: Filter type not within polyFilterCell() options.")
    }

    return(filterDF)
}

#' Polygon filtering function using "sp" package
#' @importFrom sp point.in.polygon
#' @keywords internal
pip <- function(points, pgn, points_x_column = 1, points_y_column = 2, pgn_x_column = 1, pgn_y_column = 2){
    print("F5: computing pips")
    # Points dataframe and polygon dataframe, each with only two columns corresponding to x and y values.
    if(nrow(pgn) == 0){
        ret_value <- rep(0,nrow(points))
    } else {
        pips <- sp::point.in.polygon(point.x = points[[points_x_column]],
                                     point.y = points[[points_y_column]],
                                     pol.x = pgn[[pgn_x_column]],
                                     pol.y = pgn[[pgn_y_column]])
        ret_value <- pips
        # Returning an array of 0/1 values, depending on rows fitting or not in the polygon.
    }
    return(ret_value)
}

# REVISAR - tomar en cuenta todos los filtros para decidir si una célula es filtrada por alguno o no.
calculateTruth <- function(filterDF, cell_unique_id_field = "ucid", truth_column = "filter", mode = "all"){
    print("F6.1: calculateTruth")
    # browser()
    # Descartar columnas que no me interesan para calcular la verdad
    drops <- c(cell_unique_id_field, truth_column)
    fDF <- filterDF[, !(names(filterDF) %in% drops), drop = FALSE]

    # Find filter types from the filterDF column names; names come from polyFilterCell()
    types <- sub(".*_(.*)", "\\1", names(filterDF)[!names(filterDF) %in% drops])
    types <- c("Additive" = 1, "Subtractive" = 2)[types]
    additive_types <- which(types == 1)     # identify yes-type columns
    subtractive_types <- which(types == 2)  # identify not-type columns
    
    # Tomar las columnas de filterDF que contienen los valores de verdad para cada filtro
    # y hacerles un all() por fila.
    filterDF$filter <- apply(X = fDF,
                             MARGIN = 1,        # iterate over rows
                             FUN = function(r){ # grab the row and convert it to an array
                                 r <- array(r)  # r <- array(c(T, F))
                                 r[is.na(r)] <- F  # NAs will be converted to FALSE

                                 r_at <- if (length(additive_types) == 0) T else r[additive_types]
                                 r_st <- if (length(subtractive_types) == 0) T else r[subtractive_types]

                                 if(mode == "all"){  # Subtractive overcomes Additive
                                     # all(r)        # Keep only if ALL conditions are true:
                                     all(any(r_at),  # Passes any of the aditive types
                                         all(r_st))  # Passes all of the subtractive types

                                 } else if(mode == "any"){  # Additive overcomes Subtractive
                                     # all(r)        # Keep if ANY condition is true:
                                     any(any(r_at),  # Passes any of the aditive types
                                         all(r_st))  # Passes all of the subtractive types

                                 } else {
                                     print("calculateTruth() mode not in options, defaulting to all()")
                                     all(r)
                                 }
                             })
    print("F6.2 calculateTruth")
    return(filterDF)
}

#' Add TRUTH column to cdata
#' @keywords internal
applyFilter <- function(cdataDF, filterDF, cell_unique_id_field = "ucid", truth_column = "filter"){
    print("F7.1 applyFilter")
    # Agregar a "cdata" una columna de TRUE/FALSE que refleje si cada célula pasó o no los filtros

    # Es básicamente un merge de ciertas columnas de cdataDF y cfilterDF, por "ucid"

    # Columnas que tiene el filtro en "filterDF" pero que tengo que sacar de "cdataDF" antes del merge.
    drops <- c(truth_column)

    .cdataDF <- base::merge(cdataDF[, !(names(cdataDF) %in% drops)],
                            filterDF[, c(cell_unique_id_field, truth_column)],
                            by = cell_unique_id_field)
    print("F7.2 applyFilter")
    return(.cdataDF)
}
