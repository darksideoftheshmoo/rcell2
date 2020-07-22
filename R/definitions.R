### DEFINITIONS ###

# Import operators
# `%do%` <- foreach::`%do%`
# `%>%` <- dplyr::`%>%`

#' convert sequence of numbers to ranges
#' @importFrom IRanges IRanges
#' @importFrom IRanges reduce
numbers_to_intervals <- function(numbers = c(1:10, 20:30)){
    print("F0.1")
    x <- IRanges::reduce(IRanges::IRanges(numbers, width = 1))
    x <- as.data.frame(x)[,-3]
    x <- apply(x, 1, function(x) paste(x, collapse = "-"))
    paste(x, collapse = ", ")
}

#' range expand https://www.rosettacode.org/wiki/Range_expansion#R
rangeExpand <- function(text = "1:2, 7-9", maxPos) {
    print("F0.2")
    text <- gsub("[^0-9,\\:\\-]", replacement = "", text)
    lst <- gsub("(\\d)[-\\:]", "\\1:", unlist(strsplit(text, ",")))
    numeritos <- unlist(sapply(lst, function (x) eval(parse(text=x))), use.names=FALSE)
    # numeritos <- numeritos[numeritos < maxPos]  # better to fail
    if (length(numeritos) > 0) as.numeric(numeritos) else 1:maxPos
}

#' Return brush vertices
square <- function(x1, y1, x2, y2){
    print("F0.3")
    return(list(
        x = c(x1, x1, x2, x2),
        y = c(y1, y2, y2, y1)
    ))
}

# To-do: hacer algo para poder armar filtros fuera de la app. ¿Cómo hacer filtros de una sola variable? El polígono no sirve.
#filterBox <- function(xvar, xmin, xmax, yvar, ymin, ymax)

#' Una función para procesar el facet string y que devuelva las variables presentes en names(pdata)
getFacetVars <- function(pdata, facetFormulaString = "pos ~ treatment"){
    print("F1")
    facetFormula <- eval(parse(text=facetFormulaString))
    facetVars <- all.vars(facetFormula)
    facetVars <- facetVars[facetVars %in% names(pdata)]
    return(facetVars)
}

#' Get a list of numbers from the "position" input string
getPositions <- function(input, numPos){
    print("F2")
    positions <- sort(strtoi(strsplit(input, split = ' |,|,,| ,|, ')[[1]]))
    if(length(positions) == 0) positions <- seq.int(from = 1, to = numPos, by = 1)
    return(positions)
}

#' Apply string filters to cdata
polyFilterApply <- function(polygon_df_list, cdata, truthMode = "all"){
    print("F3")
    # Initialize empty cfilter
    cfilter <- data.frame(ucid = cdata[,"ucid"], filter = T)

    # Populate cfilter columns
    for (i in seq_along(polygon_df_list)) cfilter = do.call(what = polyFilterCell,
                                                            args = list(cdataDF = cdata,
                                                                        filterDF = cfilter,  # Iterates over the output
                                                                        polygonDF = polygon_df_list[[i]],
                                                                        polygonName = paste0("polygon", i)))

    # Recalcular la verdad
    cfilter <- calculateTruth(cfilter, mode = truthMode)

    # Add TRUTH column to cdata
    cdata <- applyFilter(cdata, cfilter)

    # Return cdata and cfilter in a list.
    return(list(cdata = cdata,
                cfilter = cfilter))
}

#' Build filterDF from polygonDF and cdataDF
polyFilterCell <- function(cdataDF, filterDF, polygonDF, polygonName = 1){
    print("F4")
    # polygonDF is NULL

    xvar <- polygonDF[1, "xvar"]
    yvar <- polygonDF[1, "yvar"]

    pips <- pip(points = cdataDF[,c(xvar, yvar)], pgn = polygonDF)

    type <- polygonDF[1, "type"]  # Poolygon filter type, as chosen at the shiny app

    if(type == "Additive") {
        filterDF[, paste0(polygonName, "_", type)] <-  as.logical(pips)

    } else if(type == "Subtractive") {
        filterDF[, paste0(polygonName, "_", type)] <- !as.logical(pips)

    } else {
        print("Filter type not within polyFilterCell() options.")
    }

    return(filterDF)
}

#' Polygon filtering function using "sp" package
#' @importFrom sp point.in.polygon
pip <- function(points, pgn, points_x_column = 1, points_y_column = 2, pgn_x_column = 1, pgn_y_column = 2){
    print("F5")
    # Points dataframe and polygon dataframe, each with only two columns corresponding to x and y values.
    if(nrow(pgn) == 0){
        return(rep(0,nrow(points)))
    } else {
        pips <- sp::point.in.polygon(point.x = points[[points_x_column]],
                                     point.y = points[[points_y_column]],
                                     pol.x = pgn[[pgn_x_column]],
                                     pol.y = pgn[[pgn_y_column]])
        return(pips)
        # Returning an array of 0/1 values, depending on rows fitting or not in the polygon.
    }
}

# REVISAR - tomar en cuenta todos los filtros para decidir si una célula es filtrada por alguno o no.
calculateTruth <- function(filterDF, cell_unique_id_field = "ucid", truth_column = "filter", mode = "all"){
    print("F6")
    # browser()
    # Descartar columnas que no me interesan para calcular la verdad
    drops <- c(cell_unique_id_field, truth_column)
    fDF <- filterDF[, !(names(filterDF) %in% drops), drop = FALSE]

    # Find filter types from the filterDF column names; names come from polyFilterCell()
    types <- sub(".*_(.*)", "\\1", names(filterDF)[!names(filterDF) %in% drops])
    types <- c("Additive" = 1, "Subtractive" = 2)[types]
    additive_types <- which(types == 1)     # identify yes-type columns
    subtractive_types <- which(types == 2)  # identify not-type columns

    # Tomar las columnas de filterDF que contienen los valores de verdad para cada filtro y hacerles un all() por fila.
    filterDF$filter <- apply(X = fDF,
                             MARGIN = 1,  # by rows
                             FUN = function(r){
                                 r <- array(r)  # r <- array(c(T, F))
                                 r[is.na(r)] <- F

                                 r_at <- if (length(additive_types) == 0) T else r[additive_types]
                                 r_st <- if (length(subtractive_types) == 0) T else {r[subtractive_types]}

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

    return(filterDF)
}

#' Add TRUTH column to cdata
applyFilter <- function(cdataDF, filterDF, cell_unique_id_field = "ucid", truth_column = "filter"){
    print("F7")
    # Agregar a "cdata" una columna de TRUE/FALSE que refleje si cada célula pasó o no los filtros

    # Es básicamente un merge de ciertas columnas de cdataDF y cfilterDF, por "ucid"

    # Columnas que tiene el filtro en "filterDF" pero que tengo que sacar de "cdataDF" antes del merge.
    drops <- c(truth_column)

    cdataDF <- merge(cdataDF[, !(names(cdataDF) %in% drops)],
                     filterDF[,c(cell_unique_id_field, truth_column)],
                     by = cell_unique_id_field)
    return(cdataDF)
}




# # Build the filter string from brush and position data, and previous filter strings or selections,
# makeFilterString <- function(positions, x, y, brush, choices, selected, variableNames, negate = F){
#
#     if(length(positions) > 0){
#         positionsStr <- paste(positions, collapse = ",")
#         positionsStr <- paste("pos", "%in% c(", positionsStr, ") & ")
#     } else {
#         positionsStr <- ""
#     }
#
#     filterString <- paste(
#         positionsStr,
#         paste(x, "<=", signif(brush$xmax, 3), " & "),
#         paste(x, ">=", signif(brush$xmin, 3), " & "),
#         paste(y, "<=", signif(brush$ymax, 3), " & "),
#         paste(y, ">=", signif(brush$ymin, 3)),
#         sep = ""
#     )
#     # Negate
#     if(negate) filterString <- paste("!(", filterString, ")")
#
#     # Append filter to reactive array, this should update the UI
#     choices <- c(choices, filterString)  # all filters + new
#     # Also conserve the selected filters
#     selected <- c(selected, filterString)  # input selected + new
#
#     # Build choice names
#     choiceNames <- makeChoiceNames(choices, variableNames = variableNames)
#
#     print(choices)
#     print(choiceNames)
#
#     return(list(
#         choices = choices,
#         selected = selected,
#         choiceNames = choiceNames
#     ))
# }

# #' Funcion copada para hacer nombres bonitos a partir de un filtro tipo select
# #'
# #' @param choices A Rcell data.frame (not the object).
# #' @param variableNames Paths a la imagen de cada posición.
# #' @return Lots of stuff.
# # @examples
# # makeChoiceNames(choices, variableNames)
# #' @import foreach
# makeChoiceNames <- function(choices, variableNames){
#     # Make a nicer name for a given string filter
#     if(length(choices) > 0){
#         # If there are filters in the list, find all cdata variables used in the filter to make the filter name.
#         choiceNames <- foreach::foreach(choice=choices, .combine = c) foreach::`%do%` {
#             variableNames[sapply(variableNames, FUN = function(x){grepl(x, choice)}, USE.NAMES = F)] %>%
#                 paste(collapse = "_")
#         }
#         choiceNumbers <- paste(rep("Filter", length(choices)), seq.int(1,length(choices),1), sep = "")
#         choiceNames <- stringr::str_remove(paste(choiceNumbers, choiceNames, sep = "-"), "-$")
#     } else {
#         choiceNames <- NULL
#     }
#     return(choiceNames)
# }

# #' Funcion de filtrado para shinyCell
# #'
# #' @param cdataDF A "cdata" dataframe from Rcell or Tidycell (not the object).
# #' @param filterDF A "cfilter" dataframe with filter data. data.frame() by default.
# #' @param condition Optional: a filter as you would use in subset (quoted if string_mode = T). NA by default.
# #' @param initialDF If TRUE, supplied filterDF is ignored and a new and empty filterDF dataframe is returned. FALSE by default.
# #' @param cell_unique_id_field The cdata column wich has te primary key. "ucid" by default)
# #' @param string_mode If TRUE, the supplied condition argument is a string (a quoted filter). FALSE by default.
# #' @return Lots of stuff.
# # @examples
# # cfilter <- filterCell(cdata, initialDF = T)
# #' @import dplyr
# filterCell <- function(cdataDF, filterDF = data.frame(), condition = NA, initialDF = F, cell_unique_id_field = "ucid", string_mode = F){
#     # Esta función toma una condición (estilo subset) y agrega una columna a "cfilter" para cada fila en "cdata".
#     # Si pasa el filtro, se asigna TRUE, o bien FALSE si no pasa.
#
#     # Si se está inicializando el filtro
#     if(initialDF) filterDF <- data.frame(ucid = cdataDF[,cell_unique_id_field], filter = T)
#
#     # Build filter according to input type
#     # https://stat.ethz.ch/R-manual/R-devel/library/base/html/Logic.html
#     if(string_mode){
#         # condition <- "substitute("a.tot < 650 & a.tot > 200"
#         if(is.na(condition)) {condition_str <- "NA"} else {condition_str <- condition}
#         print(paste("Filter (string):",condition_str))
#         condition <- parse(text=condition)
#     } else {
#         # condition <- substitute(a.tot < 650 & a.tot > 200)
#         condition <- substitute(condition)
#         condition_str <- deparse(condition)
#         print(paste("Filter (expr):", condition_str))
#     }
#
#     # Si hay alguna condición asignada
#     if(condition_str != "NA"){
#         # https://stackoverflow.com/a/4605268
#         # Asignar todo true y solo poner el false. Creo que no se puede, porque no puedo asignar valores a un subset()
#
#         f <- subset(cdataDF, eval(condition))[, cell_unique_id_field, drop = FALSE]  # Tomar las que cumplen y poner TRUE
#         if(nrow(f)>0) f[, condition_str] <- TRUE
#
#         g <- subset(cdataDF, !eval(condition))[, cell_unique_id_field, drop = FALSE]  # Tomar las que NO cumplen y poner FALSE
#         if(nrow(g)>0) g[, condition_str] <- FALSE
#
#         fg <- dplyr::bind_rows(f, g)  # Juntarlas de nuevo
#         fg <- fg[order(fg$ucid),]  # y ordenarlas por "ucid"
#
#         # Mergearlas con el filterDF proporcionado en los argumentos de la función (o el inicializado al principio)
#         filterDF <- merge(filterDF,
#                           fg,
#                           by = cell_unique_id_field)
#
#         filterDF <- calculateTruth(filterDF)  # Recalcular la verdad
#     }
#
#     return(filterDF)
# }
#
# # Apply string filters to cdata
# filterCellApply <- function(stringFilters, cdata){
#     args <- stringFilters  # Test: args = c("pos %in% c( 1,2 )", "a.tot <= 446", "a.tot >= 244")
#
#     # Initialize empty cfilter
#     cfilter <- filterCell(cdata, initialDF = T)
#
#     # Populate cfilter columns
#     for (i in seq_along(args)) cfilter = do.call(what = filterCell, args = list(cdataDF = cdata,
#                                                                                 filterDF = cfilter,  # Recursive
#                                                                                 condition = args[i],
#                                                                                 string_mode = T))
#     # Add TRUTH column to cdata
#     cdata <- applyFilter(cdata, cfilter)
#
#     # Return cdata and cfilter in a list.
#     return(list(cdata = cdata, cfilter= cfilter))
# }

# MISC
# Para hacer muchos assignments en una línea
# https://strugglingthroughproblems.wordpress.com/2010/08/27/matlab-style-multiple-assignment-in%C2%A0r/

# '%=%' = function(l, r, ...) UseMethod('%=%')
#
# l = function(...) {
#     List = as.list(substitute(list(...)))[-1L]
#     class(List) = 'lbunch'
#
#     List
# }
#
# '%=%.lbunch' = function(l, r, ...) {
#     Envir = as.environment(-1)
#
#     for (II in 1:length(l)) {
#         do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
#     }
# }

# l(a,b) %=% list("2", 3)
# También está esto: https://stackoverflow.com/a/45329855/11524079

# Assign to parent dataframe
# https://stackoverflow.com/questions/42873592/assign-variable-in-parent-environment-of-a-function
#`%<-1%` <- function(x, y) { p <- parent.frame() ; p[[deparse(substitute(x))]] <- y }
#`%<-2%` <- function(x, y) { assign(deparse(substitute(x)), y, env = parent.frame())}
