#*************************************************************************#
#' cimage
#'
#' returns a cell.image object
#'
#' @param X cell.data object
#' @param ... arguments for \code{\link{cimage.cell.image}} and \code{\link{cimage.cell.data}}
#'
#' @return
#' @export
#'
#' @examples
cimage <- function(X,...) {
    UseMethod("cimage")
    }

#*************************************************************************#
#public
#ToDo: allow for transformation funcions on images
#ToDo: allow for anotation funcions, copy code from old package
#' cimage.cell.data
#'
#' calls cimage.cell.image
#'
#' @param X cell.data or cell.image object to plot
#' @param formula formula of the form 'var1+var2~var3' specifying how the images
#'  are to be ordered. See details.
#' @param facets formula of the form 'var1+var2~var3' specifying how to facet the plot. See details.
#' @param qc.filter a boolean value indicating if the quality control filter should be applied over the data
#' @param time.var variables that indicate time and should be excluded from the grouping variables.
#'  See \code{\link{get_cell_image}}
#' @param time.course boolean indicating if the image layout represents a time
#'  course and several images of the same cell at different times are expected
#' @param select character vector defining further variables that are required for the plot
#' @param exclude character vector defining variable names to be excluded
#' @param normalize.group variable names that define groups of images that
#'  should be normalized together
#' @param invert.lut boolean indicating if Look Up Table should be inverted
#' @param N Number of random cells to select from each group. If NA or 'all',
#' all cells are selected.
#' @param ... further arguments for methods. \code{cimage} calls \code{\link{get_cell_image}},
#'  so all the arguments of this function are available.
#'
#' @return
#' @export
#'
#' @examples
cimage.cell.data <- function(X,
                             formula = NULL,
                             facets = NULL,
                             qc.filter = TRUE,
                             time.var = c("*time*", "t.frame",
                                          "z.scan", "z.slice"),
                             time.course = NULL,
                             select = NULL,
                             exclude = NULL,
                             normalize.group = "channel",
                             invert.lut = FALSE,
                             N = NULL,
                             ...) {

    #defining groups
    group = c()

    if (!is.null(formula)) {
        # Rcell::.parse.formula
        pformula = .parse.formula(formula)
        group = c(pformula$lterm, pformula$rterm, pformula$sterm)
    }

    if (!is.null(facets)) {
        # Rcell::.parse.formula
        pfacets = .parse.formula(facets)
        group = c(group, pfacets$lterm, pfacets$rterm, pfacets$sterm)
    }


    #evaluating expressions in formulas
    if (isTRUE(qc.filter)) {
        X$data <- X$data[X$data$qc, ]
    }

    if (length(group) > 0) {
        for (i in 1:length(group)) {
            vars <- all.vars(parse(text = group[i]))

            if (all(vars %in% names(X$data)) & !(group[i] == vars[1])) {

                #only for expressions (to be evaluated) with all variables in X$data
                X$data[[make.names(group[i])]] <- eval(parse(text = group[i]),
                                                       X$data,
                                                       parent.frame(n = 2))
                group[i] <- make.names(group[i])
            }
        }
    }

    X$variables$all <- names(X$data)

    time.var <- base::intersect(group, .select(X$variables, time.var, warn = FALSE))

    if (isTRUE(time.var)) time.var <- NULL

    sample.cells <- any(group %in% c("...", "sample", "cell", "cells", "cellID", "ucid"))

    group <- base::setdiff(group, c(".",
                                    "...",
                                    "sample",
                                    "cell",
                                    "cells",
                                    "cellID",
                                    "ucid",
                                    "channel",
                                    time.var))

    #cell.image transformations required variables
    var_names = c()
    if (!is.null(select) | !is.null(exclude)) {

        var_names = union(var_names, .select(X$variables, select, exclude))
        var_names <- union(var_names, time.var)
        # Rcell::cnormalize
        var_names <- union(var_names, cnormalize(NULL, normalize.group))
    }

    if (length(var_names) == 0) {
        var_names = NULL
    }

    if (is.null(time.course)) {
        time.course = isTRUE(length(time.var) > 0)
    }

    #deciding how many cells to sample
    if (is.null(N)) {
        N <- 1
        if (sample.cells) {
            N <- 7
            message("sampling N=7 cells for each condition")
        }
    }

    # Rcell::get_cell_image.cell.data
    ci <- get_cell_image.cell.data(X,
                                   group = group,
                                   select = var_names,
                                   exclude = NULL,
                                   time.course = time.course,
                                   N = N,
                                   qc.filter = qc.filter,
                                   ...)

    # return(ci)

    # g: aca se rompe
    #cell.image transformations, Rcell::cnormalize
    ci <- cnormalize(ci, normalize.group)

    if (invert.lut) {
        # # Rcell::invert.lut
        ci <- invert.lut(ci)
    }

    # Rcell::cimage.cell.cimage
    return(cimage.cell.image(ci,
                             formula = formula,
                             facets = facets,
                             allow.expressions = TRUE,
                             ...))
}


#*************************************************************************#
#public
plot.Image <- function(x,
                       width = NULL,
                       height = NULL,
                       omi = 1,
                       interpolate = FALSE,
                       vp = NULL,
                       ...){

    if (is.null(width) | is.null(height)){
        ds <- dev.size()
        width = ds[1]
        height = ds[2]
    }

    img.w = dim(x)[1]
    img.h = dim(x)[2]
    w.npc = omi
    h.npc = omi

    if (img.w / img.h > width / height) {
        h.npc = w.npc*img.h/img.w*width/height
    } else {
        w.npc = h.npc * img.w / img.h * height / width
    }

    grid::grid.raster(x,
                      x = 0.5,
                      y = 0.5,
                      width = w.npc,
                      height = h.npc,
                      interpolate = interpolate,
                      vp = vp,
                      ...)
}

#*************************************************************************#
#public
if (getRversion() >= "2.15.1") utils::globalVariables(c("img.index",
                                                       "img_id_x",
                                                       "img_id_y",
                                                       "facet_x",
                                                       "facet_y",
                                                       "t.frame",
                                                       "img_xpos_start",
                                                       "facet_xpos_start",
                                                       "img_xpos_end",
                                                       "img_ypos_start",
                                                       "facet_ypos_start",
                                                       "img_ypos_end",
                                                       "channel",
                                                       "type"))

#' cimage.cell.image
#'
#' returns cropped cell images arranged in some way
#'
#' @param X cell.data or cell.image object to plot
#' @param formula formula of the form 'var1+var2~var3' specifying how the images are to be ordered. See details.
#' @param subset logical expression indicating elements or rows to keep. Don't specify channel here
#' @param facets formula of the form 'var1+var2~var3' specifying how to facet the plot. See details.
#' @param scales either 'none', 'fixed' or 'free' axis for each facet
#' @param allow.expressions allow expressions in formulas, set to TRUE when called from cimage.cell.data
#' @param nx number of columns of images within each facet. Used with \code{formula} '~var1' or 'var1~.'
#' @param ny number of rows of images within each facet. Used with \code{formulas} '~var1' or 'var1~.'
#' @param facets.nx number of columns of facets. Used with \code{facets} '~var1' or 'var1~.'
#' @param facets.ny number of rows of facets. Used with \code{facets} '~var1' or 'var1~.'
#' @param bg.col The background color of the plot
#' @param border the width in pixels of the border between images
#' @param facets.border the width in pixels of the border between facets
#' @param rev.y boolean indicating if the y axis should be reversed
#' @param font.size The size of the font to use, in pixels
#' @param font.col The color of the font to use
#' @param display boolean indicating if the created image should be displayed
#' @param ... further arguments for methods. \code{cimage} calls \code{\link{get_cell_image}},
#'  so all the arguments of this function are available.
#'
#' @return
#' @export
#'
#' @examples
cimage.cell.image <- function(X,
                              formula = NULL,
                              subset = NULL,
                              facets = NULL,
                              scales = "fixed",
                              allow.expressions = FALSE,
                              nx = NULL,
                              ny = NULL,
                              facets.nx = NULL,
                              facets.ny = NULL,
                              bg.col="white",
                              border = 1,
                              facets.border = 1,
                              rev.y = TRUE,
                              font.size = 14,
                              font.col = "black",
                              display = interactive(),
                              ...) {

    subs.subset <- substitute(subset)

    if (!is.null(subs.subset)) {
        X <- subset.cell.image(X, deparse(subs.subset))
    }

    # Rcell::img.desc
    imgdf <- droplevels(img.desc(X))

    if (any(is.na(imgdf))) {
        warning("NAs in img.desc variables: ",
                toString(names(which(sapply(imgdf, function(x) any(is.na(x)))))))
    }

    #cheking formula
    if (is.null(formula)) { #no formula
        if (is.null(facets)) { #no formula nor facets
            #just tiling images toghether
            nx = ceiling(sqrt(length(X)))
            imgdf <- transform(imgdf,
                               facet_id = 1,
                               facet_x = 1,
                               facet_y = 1,
                               facet_z = 1,
                               img_id = img.index,
                               img_id_x = 1 + (img.index - 1) %% nx,
                               img_id_y = 1 + (img.index - 1) %/% nx)

            imgdf <- transform(imgdf, img_x = img_id_x, img_y = img_id_y)

            class(X) <- "list" #this line is required to correct bug in BioC 2.13
            #g: EBImage::tile, EBImage::normalize, EBImage::combine exist
            outimg <- EBImage::tile(EBImage::normalize(EBImage::combine(X)), nx = nx)

            class(X) <- c("cell.image", "list") #this line is required to correct bug in BioC 2.13

            if (display) {
                # g: EBImage::display exists
                EBImage::display(outimg)
            }

            img.desc(outimg) <- imgdf

            return(invisible(outimg)) #like print.cell.image

        } else {#if facets but no formula, just tile images
            pformula = list(lterm = ".",
                            rterm = "img.index",
                            sterm = NULL,
                            type = "wrap_horizontal")

            if (scales != "none") {
                scales = "facets_only"
            }
        }

    } else {#formula
        if (!plyr::is.formula(formula)) {
            stop("formula required")
        }

        if (!allow.expressions) {
            # Rcell::.check.formula
            .check_formula(formula, names(imgdf))
        }

        # Rcell::.parse.formula
        pformula = .parse.formula(formula, force.valid.names = allow.expressions)

    }

    #working out facets
    if (is.null(facets) & is.null(pformula$sterm)) {#no facets

        imgdf <- transform(imgdf,
                           facet_id = factor(1),
                           facet_x = 1,
                           facet_y = 1,
                           facet_z = 1)

        pfacets <- list(lterm = ".", rterm = ".", sterm = NULL, type = "none")

        if (scales != "none") {
            scales = "fixed"
        }

    } else {#faceting

        if (!is.null(facets)) { #facet formula preset in call

            if (!plyr::is.formula(facets)) {
                stop("facets should be a formula")
            }

            if (!allow.expressions) {
                .check_formula(facets, names(imgdf))
            }

            # Rcell::.parse.formula
            pfacets <- .parse.formula(facets, force.valid.names = allow.expressions)

            pfacets$sterm <- unique(c(pformula$sterm, pfacets$sterm))

        } else {#no facet fomula, but sterm present
            pfacets <- list(lterm = ".",
                            rterm = ".",
                            sterm = pformula$sterm,
                            type = "none")
        }

        if (pfacets$type == "none" & scales != "none" & pformula$type == "grid") {
            scales = "free"
        }

        # g: Rcell::.append.panel.idxyz
        imgdf <- .append.panel.idxyz(imgdf,
                                     pfacets,
                                     var.prefix = "facet_",
                                     nx = facets.nx,
                                     ny = facets.ny,
                                     rev.y = rev.y)

    }

    # grouping vars
    group.var = setdiff(c(pformula$lterm,
                          pformula$rterm,
                          pfacets$lterm,
                          pfacets$rterm,
                          pfacets$sterm),
                        ".")

    #doing some assertions
    if (anyDuplicated(group.var)) {
        warning("some variable is duplicated between formula and facets")
    }

    if (max(plyr::ddply(imgdf, group.var, function(df) {
        data.frame(n = dim(df)[1])})$n) > 1) {
        warning("formula and facets don't uniquely specify a cell image,
                some selected images are not shown")
    }

    #calculating images positions
    if (scales %in% c("fixed", "none")) {
        # g: Rcell::.append.panel.idxyz
        imgdf <- .append.panel.idxyz(imgdf,
                                     pformula,
                                     var.prefix = "img_",
                                     nx = nx,
                                     ny = ny,
                                     rev.y = rev.y)

    } else if (scales %in% c("free",
                             "facets_only")) {
        imgdf <- do.call("rbind",
                         plyr::dlply(imgdf,
                                     plyr::.(facet_id),
                                     .append.panel.idxyz,
                                     pform = pformula,
                                     var.prefix = "img_",
                                     nx = nx,
                                     ny = ny,
                                     rev.y = rev.y))

    } else stop("unknown value for scales argument")

    #cheking all images are the same size
    img.size <- unique(sapply(X, function(x) dim(x)[1]))

    if (length(img.size) > 1) {
        stop("images of different size")
    }

    imgb <- img.size + border

    facets.df <- plyr::ddply(imgdf, plyr::.(facet_id), function(df) {
        df[1, base::intersect(names(df), c(paste("facet_",
                                                 c("id","id_x","id_y","x","y","z"),
                                                 sep = ""),
                                           pfacets$lterm,
                                           pfacets$rterm,
                                           pfacets$sterm))]})

    facets.df <- transform(facets.df, dim.x = NA, dim.y = NA)

    #building facets
    imgdf <- transform(imgdf,
                       img_xpos_start = NA,
                       img_xpos_end = NA,
                       img_ypos_start = NA,
                       img_ypos_end = NA)

    facets.img <- list()

    for (f in unique(facets.df$facet_id)) {

        df <- imgdf[imgdf$facet_id == f,]
        max.x <- max(df$img_x)
        max.y <- max(df$img_y)


        facets.df[facets.df$facet_id == f, "dim.x"] <- max.x * imgb + border
        facets.df[facets.df$facet_id == f, "dim.y"] <- max.y * imgb + border

        for (i in df$img.index) {
            x <- df[df$img.index == i, "img_x"]
            y <- df[df$img.index == i, "img_y"]

            imgdf[imgdf$facet_id == f & imgdf$img.index == i, "img_xpos_start"] <- ((x - 1) * imgb + border + 1)
            imgdf[imgdf$facet_id == f & imgdf$img.index == i, "img_xpos_end"] <- (x * imgb)
            imgdf[imgdf$facet_id == f & imgdf$img.index == i, "img_ypos_start"] <- ((y - 1) * imgb + border + 1)
            imgdf[imgdf$facet_id == f & imgdf$img.index == i, "img_ypos_end"] <- (y * imgb)
        }

    }

    #creo base de datos y lista asociada con los ejes y sus posiciones
    #modificar codigo de abajo para completar estas dos cosas
    axis.df <- data.frame()
    axis.list <- list()

    #agregar axis aca si scale free
    if (scales == "free") {
        xaxis <- plyr::dlply(imgdf, plyr::.(facet_id), .axis,
                             pform = pformula,
                             side = 1,
                             img.size = img.size,
                             border = border,
                             font.size = font.size,
                             bg.col = bg.col)#,font.col=font.col,line.col=0)

        yaxis <- plyr::dlply(imgdf, plyr::.(facet_id), .axis,
                             pform = pformula,
                             side = 2,
                             img.size = img.size,
                             border = border,
                             font.size = font.size,
                             bg.col = bg.col)#,font.col=font.col,line.col=0)

        for (i in as.character(facets.df$facet_id)) {
            xdim <- facets.df[facets.df$facet_id == i, "dim.x"] + dim(yaxis[[i]])[1]
            ydim <- facets.df[facets.df$facet_id == i, "dim.y"] + dim(xaxis[[i]])[2]

            imgdf[imgdf$facet_id == i, "img_xpos_start"] <-
                imgdf[imgdf$facet_id == i, "img_xpos_start"] + (dim(yaxis[[i]])[1] + 1)

            imgdf[imgdf$facet_id == i, "img_xpos_end"] <-
                imgdf[imgdf$facet_id == i, "img_xpos_end"] + (dim(yaxis[[i]])[1] + 1)

            axis.df <- rbind(axis.df,
                             data.frame(type = "free_xaxis",
                                        facet_id = i,
                                        facet_z = facets.df[facets.df$facet_id == i,"facet_z"],
                                        img_xpos_start = (dim(yaxis[[i]])[1] + 1),
                                        img_xpos_end = xdim,
                                        img_ypos_start = (facets.df[facets.df$facet_id == i,"dim.y"] + 1),
                                        img_ypos_end = ydim))

            axis.df <- rbind(axis.df,
                             data.frame(type = "free_yaxis",
                                        facet_id = i,
                                        facet_z = facets.df[facets.df$facet_id == i,"facet_z"],
                                        img_xpos_start = 1,img_xpos_end = (dim(yaxis[[i]])[1]),
                                        img_ypos_start = 1,img_ypos_end = (ydim - dim(xaxis[[i]])[2])))

            facets.df[facets.df$facet_id == i, "dim.x"] <- xdim
            facets.df[facets.df$facet_id == i, "dim.y"] <- ydim
        }

        names(xaxis) <- paste0("free_xaxis_", names(xaxis))
        names(yaxis) <- paste0("free_yaxis_", names(yaxis))
        axis.list <- c(axis.list, xaxis, yaxis)
    }

    #agregar facets headers aca si facet wrap
    if (pfacets$type %in% c("wrap_horizontal", "wrap_vertical")) {
        for (i in as.character(facets.df$facet_id)) {

            xdim <- facets.df[facets.df$facet_id == i, "dim.x"] #dim(facets.img[[i]])[1]
            ydim <- font.size + facets.border + facets.df[facets.df$facet_id == i, "dim.y"] #dim(facets.img[[i]])[2]
            # g: raro, no hay funcion .nchar
            xlabelpos <- (xdim / 2) - .nchar(i) * floor(font.size / 4)

            axis.df <- rbind(axis.df,
                             data.frame(type = pfacets$type,
                                        facet_id = i,
                                        facet_z = facets.df[facets.df$facet_id == i, "facet_z"],
                                        img_xpos_start = 1,
                                        img_xpos_end = xdim,
                                        img_ypos_start = 1,
                                        img_ypos_end = (font.size + facets.border)))

            facets.df[facets.df$facet_id == i, "dim.y"] <-
                facets.df[facets.df$facet_id == i, "dim.y"] + font.size + facets.border

            imgdf[imgdf$facet_id == i,"img_ypos_start"] <-
                imgdf[imgdf$facet_id == i,"img_ypos_start"] + font.size + facets.border

            imgdf[imgdf$facet_id == i, "img_ypos_end"] <-
                imgdf[imgdf$facet_id == i, "img_ypos_end"] + font.size + facets.border

            header <- EBImage::Image(bg.col,
                                     colormode = "Grayscale",
                                     dim = c(xdim, (font.size + facets.border)))

            header <- drawText(header,
                               labels = i,
                               x = xlabelpos,
                               y = font.size,
                               adj = c(0,0))

            axis.list[[paste(pfacets$type, i, sep = "_")]] <- header
        }
    }


    fct.x <- max(facets.df$dim.x) + facets.border
    fct.y <- max(facets.df$dim.y) + facets.border
    max.x <- max(facets.df$facet_x)
    max.y <- max(facets.df$facet_y)
    max.z <- max(facets.df$facet_z)

    outimg_x_px <- max.x*(fct.x + facets.border) + facets.border
    outimg_y_px <- max.y*(fct.y + facets.border) + facets.border

    #compaginar facets

    facets.df <- transform(facets.df,
                           facet_xpos_start = ((facet_x - 1) * (fct.x + facets.border) + facets.border + 1),
                           facet_ypos_start = ((facet_y - 1) * (fct.y + facets.border) + facets.border + 1))


    # g: base::subset.data.frame
    imgdf <- plyr::join(imgdf,
                  subset.data.frame(facets.df),
                  by = c("facet_id", "facet_x", "facet_y", "facet_z"))

    imgdf <- transform(imgdf,
                       img_xpos_start = img_xpos_start + facet_xpos_start,
                       img_xpos_end = img_xpos_end + facet_xpos_start,
                       img_ypos_start = img_ypos_start + facet_ypos_start,
                       img_ypos_end = img_ypos_end + facet_ypos_start)

    if (dim(axis.df)[1] > 0) {
        axis.df <- plyr::join(axis.df,
                        subset.data.frame(facets.df,
                                          select = c("facet_id",
                                                     "facet_xpos_start",
                                                     "facet_ypos_start")),
                        by = c("facet_id"))

        axis.df <- transform(axis.df,
                             img_xpos_start = img_xpos_start + facet_xpos_start,
                             img_xpos_end = img_xpos_end + facet_xpos_start,
                             img_ypos_start = img_ypos_start + facet_ypos_start,
                             img_ypos_end = img_ypos_end + facet_ypos_start)

        axis.df$facet_xpos_start <- NULL
        axis.df$facet_ypos_start <- NULL
    }


    #agregar facets axis si facet grid
    h.xfacetsaxis = 0
    if (pfacets$type %in% c("grid") & scales != "none") {

        facets.df <- transform(facets.df, img_x = facet_x, img_y = facet_y)

        # g: Rcell::.axis
        xfacetsaxis <- .axis(facets.df,
                             pform = pfacets,
                             side = 3,
                             img.size = fct.x,
                             border = facets.border,
                             font.size = font.size,
                             bg.col = bg.col)#,font.col=font.col,line.col=0)

        yfacetsaxis <- .axis(facets.df,
                             pform = pfacets,
                             side = 4,
                             img.size = fct.y,
                             border = facets.border,
                             font.size = font.size,
                             bg.col = bg.col)#,font.col=font.col,line.col=0)

        h.xfacetsaxis <- dim(xfacetsaxis)[2]

        imgdf <- transform(imgdf,
                           img_ypos_start = img_ypos_start + h.xfacetsaxis,
                           img_ypos_end = img_ypos_end + h.xfacetsaxis)

        if (dim(axis.df)[1] > 0) {
            axis.df <- transform(axis.df,
                                 img_ypos_start = img_ypos_start + h.xfacetsaxis,
                                 img_ypos_end = img_ypos_end + h.xfacetsaxis)
        }

        for (z in 1:max.z) {
            axis.df <- rbind(axis.df,
                             data.frame(type = "grid_xfacetaxis",
                                        facet_id = "all",
                                        facet_z = z,
                                        img_xpos_start = 1,
                                        img_xpos_end = dim(xfacetsaxis)[1],
                                        img_ypos_start = 1,
                                        img_ypos_end = h.xfacetsaxis))

            axis.df <- rbind(axis.df,
                             data.frame(type = "grid_yfacetaxis",
                                        facet_id = "all",
                                        facet_z = z,
                                        img_xpos_start = (outimg_x_px + 1),
                                        img_xpos_end = (outimg_x_px + dim(yfacetsaxis)[1]),
                                        img_ypos_start = (h.xfacetsaxis + 1),
                                        img_ypos_end = (outimg_y_px + h.xfacetsaxis)))
        }

        axis.list[["grid_xfacetaxis_all"]] <- xfacetsaxis
        axis.list[["grid_yfacetaxis_all"]] <- yfacetsaxis

        outimg_x_px <- outimg_x_px + dim(yfacetsaxis)[1]
        outimg_y_px <- outimg_y_px + dim(xfacetsaxis)[2]
    }

    #agregar axes aca si scale fixed y formula grid
    if (scales == "fixed" & pformula$type == "grid") {

        xaxis <- .axis(imgdf,
                       pform = pformula,
                       side = 1,
                       img.size = img.size,
                       border = border,
                       font.size = font.size,
                       bg.col = bg.col)#,font.col=font.col,line.col=0)

        yaxis <- .axis(imgdf,
                       pform = pformula,
                       side = 2,
                       img.size = img.size,
                       border = border,
                       font.size = font.size,
                       bg.col = bg.col)#,font.col=font.col,line.col=0)


        axis.list[["fixed_xaxis_all"]] <- xaxis
        axis.list[["fixed_yaxis_all"]] <- yaxis

        xdim <- outimg_x_px + dim(yaxis)[1]
        ydim <- outimg_y_px + dim(xaxis)[2]

        imgdf <- transform(imgdf,
                           img_xpos_start = img_xpos_start + (dim(yaxis)[1]),
                           img_xpos_end = img_xpos_end + (dim(yaxis)[1]))

        if (dim(axis.df)[1] > 0) {
            axis.df <- transform(axis.df,
                                 img_xpos_start = img_xpos_start + (dim(yaxis)[1]),
                                 img_xpos_end = img_xpos_end + (dim(yaxis)[1]))
        }

        # tmp<-list()
        for (z in 1:max.z) {

            for (i in 1:max.y) {
                axis.df <- rbind(axis.df,
                                 data.frame(type = "fixed_yaxis",
                                            facet_id = "all",
                                            facet_z = z,
                                            img_xpos_start = 1,
                                            img_xpos_end = (dim(yaxis)[1]),
                                            img_ypos_start = (i * (fct.y + facets.border) - dim(yaxis)[2] + 1 + h.xfacetsaxis),
                                            img_ypos_end = (i * (fct.y + facets.border) + h.xfacetsaxis)))


            }

            for (i in 1:max.x) {
                axis.df <- rbind(axis.df,
                                 data.frame(type = "fixed_xaxis",
                                            facet_id = "all",
                                            facet_z = z,
                                            img_xpos_start = (dim(yaxis)[1] + 1 + (i - 1) * (fct.x) + i * facets.border),
                                            img_xpos_end = (dim(yaxis)[1] + (i - 1) * (fct.x) + i * facets.border + dim(xaxis)[1]),
                                            img_ypos_start = (outimg_y_px + 1), img_ypos_end = ydim))


            }
        }

        outimg_x_px <- outimg_x_px + dim(yaxis)[1]
        outimg_y_px <- outimg_y_px + dim(xaxis)[2]

    }

    # g: EBImage::Image exists
    outimg2 <- plyr::rlply(max.z,
                           EBImage::Image(bg.col,
                                          colormode = "Grayscale",
                                          dim = c(outimg_x_px,outimg_y_px)))

    for (i in imgdf$img.index) {
        df <- imgdf[imgdf$img.index == i,]
        outimg2[[df$facet_z]][df$img_xpos_start:df$img_xpos_end, df$img_ypos_start:df$img_ypos_end] <- X[[i]]
    }

    if (dim(axis.df)[1] > 0) {
        axis.df <- transform(axis.df, img_id = paste0(type, "_", facet_id))
        for (i in 1:dim(axis.df)[1]) {
            df <- axis.df[i,]
            outimg2[[df$facet_z]][df$img_xpos_start:df$img_xpos_end, df$img_ypos_start:df$img_ypos_end] <- axis.list[[as.character(df$img_id)]]
        }
    }

    # g: EBImage::combine exists
    outimg2 <- EBImage::combine(outimg2)


    imgdf <- plyr::rename(imgdf,
                          c("img_xpos_start" = "xpos_start",
                            "img_xpos_end" = "xpos_end",
                            "img_ypos_start" = "ypos_start",
                            "img_ypos_end" = "ypos_end",
                            "facet_z" = "slice"))

    if ("t.frame" %in% names(imgdf)) {
        imgdf <- plyr::arrange(imgdf, pos, cellID, channel, t.frame)
    } else {
        imgdf <- plyr::arrange(imgdf, pos, cellID, channel)
    }

    attr(outimg2, "img.desc") <- imgdf

    if (display) {
        EBImage::display(outimg2, title = "cimage")
    }

    EBImage::display(outimg2, title = "cimage")
    #return(invisible(outimg2))
    }

#*************************************************************************#
#public
cimage.default <- function(X, ...){
    X <- get_cell_image.default(X, ...)
    return(cimage.cell.image(X, ...))
}

#*************************************************************************#
#generic
#returns a cell.image object
#' Title
#'
#' @param X cell.data object or data.frame that specifies the images
#' @param ... further arguments for methods
#'
#' @return
#' @export
#'
#' @examples
get_cell_image <- function(X,...) {
    UseMethod("get_cell_image")
    }

#*************************************************************************#
#public
#returns a list of croped cells images, given a cell.data object
#
#' get_cell_image
#'
#'
#' Retrieves the images from single cells in an cell.image object
#'
#' @param X cell.data object or data.frame that specifies the images
#' @param subset logical expression indicating elements or rows to keep. Don't specify channel here.
#' @param channel.subset logical expression to specify which image to retrieve with channel and t.frame variables.
#' @param channel character vector of channels to retrieve. If specified, defines the order of the channels.
#' @param time.course boolean indicating if the desired image montage is a time course (i.e. several images for the same cell)
#' @param group character vector or quoted names of variables who's interaction define the groups from which select \code{N} random cells.
#' @param na.rm boolean indicating if NAs should be removed.
#' @param N Number of random cells to select from each group. If NULL all cells are selected
#' @param select character vector defining variables names to be included in the returned cell.image object
#' @param exclude character vector defining variables names to be excluded from the returned cell.image object
#' @param qc.filter a boolean value indicating if the quality control filter should be applied over the data
#' @param box.size size in pixels of the image containing the cells. This specifies the 'radius', i.e. the image will be a square of length 2*box.size+1
#' @param ... further arguments for methods
#'
#' @return
#' @export
#'
#' @examples
get_cell_image.cell.data <- function(X,
                                     subset = NULL,
                                     channel.subset = NULL,
                                     channel = NULL,
                                     time.course = TRUE,
                                     group = NULL,
                                     na.rm = TRUE,
                                     N = 7,
                                     select = NULL,
                                     exclude = NULL,
                                     qc.filter = TRUE,
                                     box.size = 20,
                                     ...){

    subset = substitute(subset)

    channel.subset = substitute(channel.subset)
    group = as.character(group)

    #filtering by qc variable
    if (class(X$data$qc) == "logical" && qc.filter) {
        data = subset(X$data,qc)
    } else {
        data = X$data
    }

    #checking that group variables exists
    if (!all(group %in% X$variables$all)) {
        stop(paste(group[!group %in% X$variables$all], collapse = ", "),
             " not found in cell.data")
    }

    #selecting variables to include in dataset
    var_names = union(all.vars(channel.subset), group)
    var_names = union(var_names, all.vars(subset))
    if ("maskID" %in% names(X$data)) {
        var_names = union(var_names,"maskID")
    }

    var_names = union(var_names, c(.CELLID_ID_VARS, "ucid", "xpos", "ypos"))

    if (!is.null(select) | !is.null(exclude)) {
        var_names = union(var_names,.select(X$variables,select,exclude))
    }
    var_names = base::setdiff(var_names, "channel")

    #subsetting the dataset
    if (!is.null(subset)) {
        data = data[eval(subset, data), base::intersect(var_names, names(data))]

    } else {
        data = subset(data, select = var_names)
    }

    #cheking for NAs in grouping vars
    if (any(is.na(data[,group]))) {
        NA.vars <- group
        if (length(group) > 1) {
            NA.vars <- names(which(sapply(data[,group], function(x) any(is.na(x)))))}

        if (length(NA.vars) > 1) {
            NA.registers <- apply(is.na(data[,NA.vars]), 1, any)
        } else {
            NA.registers <- is.na(data[,NA.vars])
        }

        if (isTRUE(na.rm)) { #removing registers with NAs in grouping vars
            data <- data[!NA.registers, ]
            message("Removing ", sum(NA.registers), " registers (",
                    signif(100 * sum(NA.registers) / length(NA.registers), 2),
                    "%) with NAs in group variable(s): ", toString(NA.vars))
        } else {
            warning("There are ", sum(NA.registers), " registers (",
                    signif(100 * sum(NA.registers) / length(NA.registers), 2),
                    "%) with NAs in group variable(s): ", toString(NA.vars),
                    ".\nUse na.rm=TRUE to remove them.")
        }
    }

    #selecting N random cells per group
    if (is.null(N) | is.na(N) | N == 0 | tolower(N) == "all") {
        N <- NULL
    }

    if (length(group) > 0) { #sampling by group
        data <- do.call("rbind", plyr::dlply(data, group, .sample.N.ucid, N))
    } else {
        #random sampling if no groups are specifyied
        data <- .sample.N.ucid(data, N)
    }

    # if several times per cell in a group AND NOT set as time.course
    # AND several t.frames from which to choose
    if (!time.course & max(table(data$ucid)) > 0 & length(unique(data$t.frame)) > 1) {
        #selecting random times
        data <- plyr::ddply(data,
                      c(group, "ucid"),
                      function(df) df[sample(1:dim(df)[1], 1), ])
        message("Selecting random t.frames for each cell. Use time.course=TRUE if it's a time course.")
    }

    #"unfolding" dataset
    X$images <- transform(X$images, image = factor(image))

    #removing vars that might cause problems when merging
    merge.vars <- base::setdiff(names(X$images), base::setdiff(names(data), c("pos", "t.frame")))

    data <- merge(data, X$images[,merge.vars], by = c("pos", "t.frame"), suffixes = c("",""))

    #filtering by channel
    if (!is.null(channel.subset)) {
        data = subset(data, eval(channel.subset, data, parent.frame(n = 1)))
    }

    if (!is.null(channel)) {
        data = data[data$channel %in% channel, ]
        #redifining channel as ordered factor
        data$channel <- factor(data$channel, levels = channel, ordered = TRUE)
    }

    return(get_cell_image.data.frame(data, box.size = box.size, ...))
}

#*************************************************************************#
#public
#returns a list of croped cells images, given a data.frame specifying
#the xpos, ypos, path and image name
if (getRversion() >= "2.15.1") utils::globalVariables(c("qc",
                                                        "path.image",
                                                        "xpos",
                                                        "ypos",
                                                        "offset.x",
                                                        "offset.y"))

#' get_cell_image.data.frame
#'
#' @param X cell.data object or data.frame that specifies the images
#' @param box.size size in pixels of the image containing the cells. This specifies the 'radius', i.e. the image will be a square of length 2*box.size+1
#' @param contained.box boolean indicating if the XY position of the box should be corrected to be contained in the original image. Relevant for cells close to the image border. If FALSE the part of the box outside the original image will be filled with \code{bg.col}
#' @param bg.col color to be used for the background of the images
#' @param ... further arguments for methods
#'
#' @return
#' @export
#'
#' @examples
get_cell_image.data.frame <- function(X,
                                      box.size = 20,
                                      contained.box = FALSE,
                                      bg.col = 0,
                                      ...){

    #renaming data.frame X to img.desc for clarity
    img.desc <- X

    #cheking for required data.frame columns
    if (!all(c("xpos", "ypos", "image", "path") %in% names(img.desc))) {
        stop("xpos, ypos, image and path columns required in data.frame")
    }

    if (!"offset.x" %in% names(img.desc)) {
        img.desc <- transform(img.desc,offset.x = 0)
    }

    if (!"offset.y" %in% names(img.desc)) {
        img.desc <- transform(img.desc,offset.y = 0)
    }

    #cheking that the image files exist
    img.fnames <- with(img.desc, paste(path, image, sep = "/"))
    img.fnames.exist <- file.exists(img.fnames)

    if (!all(img.fnames.exist)) {
        stop("image file ", img.fnames[!img.fnames.exist][1], " not found")
    }

    #adding image index and center.shift vars
    img.desc <- transform(img.desc,
                          img.index = seq(1, length(img.desc$image)),
                          path.image = paste(img.desc$path, img.desc$image, sep = "/"),
                          center.x.shift = NA,
                          center.y.shift = NA)

    #warning if it is going to take long
    img.to.open <- length(unique(img.desc$path.image))
    cell.to.crop <- length(unique(img.desc$img.index))

    if (img.to.open > 100 | cell.to.crop > 1000) {
        message(img.to.open,
                " images will be opened and ",
                cell.to.crop,
                " cells will be cropped. This can take some time.")
    }

    # open each image once, and crop all the required cells, then close image
    cell.image <- list()
    for (i in as.character(unique(img.desc$path.image))) {
        df <- subset(img.desc, path.image == i)
        is.mask <- isTRUE(unique(df$channel) == "mask")
        if (is.mask) { #mask are saved as .rds objects
            img <- readRDS(i)
        } else {#regular tif image
            #capture.output to avoid anoying text
            # g: EBImage::readImage holi_2019-03-05
            capture.output(img <- EBImage::readImage(i))

            if (EBImage::colorMode(img) == EBImage::Color) {
                EBImage::colorMode(img) <- EBImage::Grayscale
                img <- img[,,1]
            }
        }

        max.x <- dim(img)[1]
        max.y <- dim(img)[2]

        #calculating cell box edges
        df <- transform(df,
                        x0 = xpos + offset.x - box.size,
                        x1 = xpos + offset.x + box.size,
                        y0 = ypos + offset.y - box.size,
                        y1 = ypos + offset.y + box.size)

        for (j in unique(df$img.index)) {

            #calculating cell box edges
            x0 <- df[df$img.index == j, "x0"]
            x1 <- df[df$img.index == j, "x1"]
            y0 <- df[df$img.index == j, "y0"]
            y1 <- df[df$img.index == j, "y1"]

            left.margin <- 0
            right.margin <- 0
            top.margin <- 0
            bottom.margin <- 0
            center.x.shift <- 0
            center.y.shift <- 0

            if (contained.box) {
                #correcting for image border

                if (x0 < 1) {
                    center.x.shift <- x0 - 1
                    x0 <- 1
                    x1 <- 2*box.size + 1
                }

                if (x1 > max.x) {
                    center.x.shift <- x1-max.x
                    x0 <- max.x - 2 * box.size
                    x1 <- max.x
                }

                if (y0 < 1){
                    center.y.shift <- y0 - 1
                    y0 <- 1
                    y1 <- 2 * box.size + 1
                }

                if (y1 > max.y) {
                    center.y.shift <- y1 - max.y
                    y0 <- max.y - 2 * box.size
                    y1 <- max.y
                }

            } else {
                if (x0 < 1) {
                    left.margin = 1 - x0
                    x0 <- 1
                }
                if (x1 > max.x) {
                    right.margin = x1 - max.x
                    x1 <- max.x
                }

                if (y0 < 1) {
                    top.margin = 1 - y0
                    y0 <- 1
                }

                if (y1 > max.y) {
                    bottom.margin = y1 - max.y
                    y1 <- max.y
                }
            }

            tmp <- EBImage::Image(bg.col,
                                  colormode = "Grayscale",
                                  dim = c(2 * box.size + 1, 2 * box.size + 1))

            tmp[(left.margin+1):(2*box.size+1-right.margin),(top.margin+1):(2*box.size+1-bottom.margin)]<-img[x0:x1,y0:y1]

            #removing masks of other cells
            if (is.mask) {
                mid <- df$maskID[df$img.index == j]
                tmp[tmp != mid] <- 0
                tmp[tmp == mid] <- 1
            }

            cell.image[[j]] <- tmp

            img.desc[img.desc$img.index == j,c("center.x.shift","center.y.shift")] <- c(center.x.shift,center.y.shift)
        }
    }

    class(cell.image) <- c("cell.image", "list")
    attr(cell.image, "img.desc") <- subset(img.desc, select = -path.image)

    return(cell.image)
}

#*************************************************************************#
# calculate the displacement of the cell center from the center of the tile genereted by contained.box=TRUE



#*************************************************************************#
#generic
#returns a cell.image object
get_cell_image.default <-  function(X,
                                    box.size = 20,
                                    ...) {

    get_cell_image.data.frame(as.data.frame(X),
                              box.size = box.size,
                              ...)
}

#*************************************************************************#
#public
#checks if an objects is a cell.data object
is.cell.image <- function(X) {
    inherits(X, "cell.image")
}

#*************************************************************************#
#public
merge.cell.image <- function(x, y, ...) {

    if (!is.cell.image(x) | !is.cell.image(y)) {
        stop("x and y should be of class cell.image")
    }
    box.size.x <- unique(as.vector(sapply(x, dim)))

    if (length(box.size.x) > 1) {
        stop("box.size for x images is not homogeneous")
    }

    box.size.y <- unique(as.vector(sapply(y, dim)))

    if (length(box.size.y) > 1) {
        stop("box.size for y images is not homogeneous")
    }

    if (box.size.x != box.size.y) {
        stop("x and y should have the same box.size")
    }

    img.desc.x <- img.desc(x)
    img.desc.y <- img.desc(y)

    if (".merge.index" %in% names(img.desc.x)) { #x has been merged before
        mi <- max(img.desc.x$.merge.index)
    } else {#first time x is merged
        img.desc.x$.merge.index <- mi <- 1
    }

    if (".merge.index" %in% names(img.desc.y)) { #y has been merged before
        img.desc.y$.merge.index <- img.desc.y$.merge.index + mi
    } else {#first time y is merged
        img.desc.y$.merge.index <- mi + 1
    }

    idn <- base::intersect(names(img.desc.x), names(img.desc.y))
    excluded.vars <- base::setdiff(union(names(img.desc.x), names(img.desc.y)), idn)

    if (length(excluded.vars) > 0) {
        warning("variables ", toString(excluded.vars), " are only in one object and will be excluded")
    }

    img.desc.x <- subset(img.desc.x, select = idn)
    img.desc.y <- subset(img.desc.y, select = idn)

    img.desc.y$img.index <- img.desc.y$img.index + length(x)
    z <- c(x, y)
    attr(z, "img.desc") <- rbind(img.desc.x, img.desc.y)
    class(z) <- c("cell.image", "list")
    return(z)
}


#*************************************************************************#
#public
#returns img.desc attribute of cell.image object
img.desc <- function(X) {
    if (!(is.cell.image(X) | EBImage::is.Image(X))) {
        stop("cell.image or Image object required")
    }
    attr(X, "img.desc")
}

#*************************************************************************#
#public
# g: replacement function! lets you replace the img.desc attribute of a
# cell.image object like this: img.desc(object) <- "whatever"
"img.desc<-"  <- function(X, value){
    attr(X, "img.desc") <- value
    X
}

#*************************************************************************#
#public
summary.cell.image <- function(object, ...) {
    summary = attr(object, "img.desc")
    class(summary) <- c("summary.cell.image", "data.frame")
    return(summary)
}

#*************************************************************************#
#public
print.summary.cell.image <- function(x, ...) {
    cat("cell images from", toString(unique(x$path)), "\n")
    cat("\npositions:", rcell2::.format.sequence(unique(x$pos)))
    cat("\ntime frames:", rcell2::.format.sequence(unique(x$t.frame)))
    cat("\nchannels: ", toString(unique(x$channel)), sep = "")
    cat("\nimage index: ", rcell2::.format.sequence(unique(x$img.index)), sep = "")
    cat("\nnumber of cells: ", length(unique(interaction(x$pos,x$cellID,drop = TRUE))), sep = "")
    cat("\n\n")

    var_names = c("img.index", "pos", "cellID", "t.frame", "channel")

    var_names = c(var_names,
                  base::setdiff(names(x), c(var_names, "image", "path", "xpos","ypos")))

    print.data.frame(head(subset(x, select = var_names)),
                     row.names = FALSE)
}


#*************************************************************************#
#public
print.cell.image <- function(x, nx = ceiling(sqrt(length(x))), ...) {
    xtitle <- paste("cell.image from", toString(unique(img.desc(x)$path)))
    class(x) <- "list"
    EBImage::display(EBImage::tile(EBImage::normalize(EBImage::combine(x)), nx = nx),
                     title = xtitle)
}

#*************************************************************************#
#public
if (getRversion() >= "2.15.1") utils::globalVariables(c("img.index"))
subset.cell.image <- function(x,
                              .subset = TRUE,
                              select = NULL,
                              exclude = NULL,
                              ...) {

    subs.subset = substitute(.subset)

    img.desc.x <- img.desc(x)
    img.desc.x <- plyr::arrange(img.desc.x, img.index)
    select.vars <- .select(list(all = names(img.desc.x)), select, exclude)

    if (isTRUE(select.vars)){
        select.vars = names(img.desc.x)}

    if (subs.subset[[1]] == "deparse") { #called from other function
        keep.rows <- eval(parse(text = .subset), img.desc(x))
    } else {
        keep.rows <- eval(subs.subset, img.desc.x)
    }

    x <- x[img.desc.x[keep.rows, "img.index"]]
    img.desc.x <- img.desc.x[keep.rows, select.vars]
    img.desc.x$img.index <- seq_len(length(x))

    class(x) <- c("cell.image", "list")
    img.desc(x) <- img.desc.x
    return(x)
}


#*************************************************************************#
#public
#ToDo: implement annotate
show_img <- function(X,
                     pos,
                     t.frame = 0,
                     channel = "BF.out",
                     image.title = "",
                     annotate = NULL,
                     cross = !qc,
                     qc.filter = FALSE,
                     subset = TRUE,
                     cross.col = c(0.1,0.9),
                     display = interactive(),
                     normalize = TRUE,
                     ...) {

    cross = substitute(cross)
    subset = substitute(subset)

    if (any(!is.element(pos,X$images$pos))) {
        stop("Selected positions not in dataset")
    }

    arg.df <- data.frame(pos = rep(pos,each = length(t.frame) * length(channel)),
                         t.frame = rep(rep(t.frame, each = length(channel)), length(pos)),
                         channel = rep(channel,length(pos) * length(t.frame)))

    # g: original: join
    img.df <- plyr::join(arg.df, X$images, by = c("pos", "t.frame", "channel"))
    img.df <- data.frame(img.df, index = 1:dim(img.df)[1])

    #cheking that the image files exist
    img.fnames = with(img.df, paste(path, image, sep = "/"))
    img.fnames.exist = file.exists(img.fnames)
    if (!all(img.fnames.exist)) {
        stop("image file ", img.fnames[!img.fnames.exist][1], " not found")
    }

    #loading the images
    img.list <- list()

    for (i in 1:length(img.fnames)) {

        #capture.output to avoid anoying text
        # g: utils::capture.output
        capture.output(img.list[[i]] <- EBImage::readImage(img.fnames[i]))

        if (isTRUE(normalize)) {
            img.list[[i]] <- EBImage::normalize(img.list[[i]])
        }
    }

    #subsetting the data
    X$data <- X$data[X$data$pos %in% pos & X$data$t.frame %in% t.frame,]

    if (dim(X$data)[1] > 0) {
        if (isTRUE(qc.filter) && class(X$data$qc) == "logical") {
            X$data = subset(X$data, qc)
        }

        if (!isTRUE(subset)) {
            X$data <- subset(X$data, eval(subset, X$data))
        }

        #adding crosses to image
        if (!is.null(cross)) {
            cross.df = data.frame(subset(X$data, select = c(pos,t.frame,cellID,xpos,ypos)),
                                  cross = as.logical(eval(cross,X$data)))

            for (i in img.df$index) {
                cell.df = cross.df[cross.df$pos == img.df$pos[img.df$index == i] &
                                       cross.df$t.frame == img.df$t.frame[img.df$index == i] &
                                       cross.df$cross,
                                   c("xpos", "ypos")]

                for (j in 1:length(cross.col)) {

                    # RCell::drawCross
                    img.list[[i]] <- drawCross(img.list[[i]],
                                               cell.df$xpos + (j - 1),
                                               cell.df$ypos,
                                               col = cross.col[j])
                }
                # display(img.list[[i]])
            }
        }

        if (!is.null(annotate)) {
            stop("annotate not implemented yet, sorry")
        }

    } else {
        if (!is.null(cross)) {
            message("no cells for this image")
        }
    }

    SHOW_IMAGE <- EBImage::combine(img.list)

    if (display) {
        EBImage::display(SHOW_IMAGE, title = "show_image")
    }

    return(invisible(SHOW_IMAGE))
}

# creates alias
show_image <- show_img

#*************************************************************************#
#public
update_img_path <- function(X,
                            img.path = getwd(),
                            subset = NULL) {

    if (!isTRUE(is.cell.data(X))) {
        stop("First argument should be of class cell.data")
    }

    #cheking that the image files exist in the new path
    img.fnames = with(X$images, paste(img.path, image, sep = "/"))
    img.fnames.exist = file.exists(img.fnames)

    if (sum(img.fnames.exist) < length(img.fnames.exist)){
        warning("some images could not be found in new path, e.g. ",
                X$images$image[!img.fnames.exist][1])
    }

    subset = substitute(subset)

    if (is.null(subset)) {
        X$images$path <- as.factor(img.path)
    } else {
        X$images$path[eval(subset, X$images, parent.frame(n = 1))] <- as.factor(img.path)
    }

    return(X)
}

#*************************************************************************#
#public
write.cell.image <- function(x, file, ...) {
    db <- img.desc(x)

    if (class(x)[1] == "cell.image") {
        class(x) <- "list"
        x <- EBImage::combine(x)
    }

    EBImage::writeImage(x, file, ...)

    if (!is.null(db)) {
        file <- paste0(sub("[.][^.]*$", "", file, perl = TRUE), "-imgdesc.txt")

        rcell2::write.delim(db, file)
    }

    return(invisible(NULL))
}

#*************************************************************************#
#public
read.cell.image <- function(file, ...) {

    file.db <- paste0(sub("[.][^.]*$", "", file, perl = TRUE),"-imgdesc.txt")

    if (!file.exists(file)) stop("file not found")

    x <- EBImage::readImage(file,...)

    if (!file.exists(file.db)) warning("image descripction file not found")

    db <- read.delim(file.db)

    if ("img.index" %in% names(db)) { #dealing with a cell.image
        x <- lapply(seq_len(EBImage::numberOfFrames(x)),
                    FUN = function(i) EBImage::getFrame(x, i))

        class(x) <- c("cell.image", "list")
    }

    attr(x, "img.desc") <- db

    return(x)
}

#*************************************************************************#
#public
revFactor <- function(x) {
    factor(x, levels = rev(levels(x)), ordered = TRUE)
}

###########################################################################
#####################cell.image transformation functions###################
###########################################################################

cnormalize <- function(X = NULL,
                       normalize.group = c("channel"),
                       ft = c(0,1),
                       ...) {

    normalize.group = as.character(normalize.group)

    if (is.null(X)) {
        return(base::setdiff(normalize.group, c("channel", "sample", "...")))
    }

    if (length(normalize.group) == 0) {
        return(X)
    }



    img.list <- plyr::dlply(img.desc(X), normalize.group, function(df) df$img.index)

    for (i in names(img.list)) {

        img <- EBImage::combine(X[img.list[[i]]])
        img <- EBImage::normalize(img, separate = FALSE, ft = ft)
        n.frames <- length(img.list[[i]])

        if (n.frames == 1) {
            X[[img.list[[i]]]] <- img
        } else {
            for (j in 1:n.frames)
                X[[img.list[[i]][j]]] <- img[,,j]
        }
    }

    return(X)
}


invert.lut <- function(X = NULL, ...){
    if (is.null(X)) {
        return(c("ucid", "channel", "t.frame"))
    }

    for (i in seq_along(X)) {
        X[[i]] <- (1 - X[[i]])
    }

    return(X)
}

#*************************************************************************#
#cell.image transformation functions
#cell.image cell.image apply, extension of plyr package for cell image
#ToDo: allow for functionss that reuturn image stacks o same length as input
#	   set MARGIN=NULL for functions like normalize and MARGIN=3 in a per image function
#	   manage img.desc appropiately


ciciply <- function(X = NULL,
                    group = c("pos", "cellID", "channel"),
                    FUN = sum,
                    MARGIN = c(1, 2),
                    warn = TRUE) {

    if (is.null(X)) {
        return(group)
    }

    img.list <- plyr::dlply(img.desc(X), group, function(df) df$img.index)
    Y <- list()
    Y.index <- data.frame(img.index = seq_len(length(img.list)), img.name = names(img.list))

    for (i in names(img.list)) {
        img <- EBImage::combine(X[img.list[[i]]])
        img <- EBImage::as.Image(apply(img, MARGIN, FUN))
        Y[[Y.index[Y.index$img.name == i, "img.index"]]] <- img
    }

    Y.img.desc <- plyr::ddply(img.desc(X), group, plyr::colwise(.unique.na))

    rm.vars <- names(Y.img.desc)[sapply(Y.img.desc, FUN = function(x) any(is.na(x)))]

    if (warn) {
        message("removing variables from img.desc: ", toString(rm.vars))
    }

    Y.img.desc <- subset(Y.img.desc, select = base::setdiff(names(Y.img.desc), rm.vars))
    Y.img.desc$img.name <- interaction(Y.img.desc[,group], drop = TRUE, lex.order = FALSE)

    if ("img.index" %in% names(Y.img.desc)) {
        Y.img.desc$img.index <- NULL
    }

    Y.img.desc <- plyr::join(Y.img.desc, Y.index, by = "img.name")
    Y.img.desc$img.name <- NULL

    class(Y) <- c("cell.image", "list")
    attr(Y, "img.desc") <- Y.img.desc

    return(Y)
}

add.nucleus.boundary <- function(X = NULL,
                                 radii = c(2,3,4,5,6,7),
                                 pos.nucl.channel = "YFP",
                                 col = 0.75,
                                 ...) {

    if (is.null(X)) {
        return(c("xpos","ypos","xpos.nucl.?","xpos.nucl.?"))
    }

    db = img.desc(X)
    xpos.nucl.var = paste("xpos.nucl.", tolower(substr(pos.nucl.channel, 1, 1)), sep = "")
    ypos.nucl.var = paste("ypos.nucl.", tolower(substr(pos.nucl.channel, 1, 1)), sep = "")

    for (i in 1:length(X)){
        img <- X[[i]]
        xcenter = ceiling(dim(X[[i]])[2] / 2) + db[i, xpos.nucl.var] - db[i, "xpos"]
        ycenter = ceiling(dim(X[[i]])[1] / 2) + db[i, ypos.nucl.var] - db[i, "ypos"]

        for (r in radii){
            img <- EBImage::drawCircle(img, xcenter, ycenter, r, col = col)
        }
        X[[i]] <- img
    }

    return(X)
}

add.maj.min.axis <- function(X = NULL,
                             col = 0.75,
                             angle.var = NA,
                             ...) {

    if (is.null(X)) {
        return(c("xpos","ypos","maj.axis","min.axis"))
    }

    db = img.desc(X)

    for(i in 1:length(X)){
        img <- X[[i]]
        xcenter = ceiling(dim(X[[i]])[2] / 2)
        ycenter = ceiling(dim(X[[i]])[1] / 2)
        angle = 0

        if(!is.na(angle.var)) {
            angle = db[i, angle.var]
        }

        majAxis = db[i, "maj.axis"] / 2
        minAxis = db[i, "min.axis"] / 2

        img <- drawLine(img,
                        round(xcenter - majAxis * cos(angle)),
                        round(ycenter + majAxis * sin(angle)),
                        round(xcenter + majAxis * cos(angle)),
                        round(ycenter - majAxis * sin(angle)))

        img <- drawLine(img,
                        round(xcenter - minAxis * sin(angle)),
                        round(ycenter - minAxis * cos(angle)),
                        round(xcenter + minAxis * sin(angle)),
                        round(ycenter + minAxis * cos(angle)))

        X[[i]] <- img
    }

    return(X)
}

#####################General Image Manipulation Functions###################

drawCross<-function(img, x, y, size=2, col=0.75, z=1){
    #EBImage::validImage(img)
    if (EBImage::colorMode(img) == EBImage::Color)
        stop("this method doesn't support the 'Color' color mode")
    if (any(is.na(img)))
        stop("'img' shouldn't contain any NAs")
    if (missing(x)) stop("'x' is missing")
    if (missing(y)) stop("'y' is missing")
    if (z < 1 | z > EBImage::numberOfFrames(img, "render"))
        stop("'z' must be a positive integer lower than the number of image frames")
    if (EBImage::colorMode(img) == EBImage::Color) {
        rgb = as.numeric(col2rgb(col)/255)
        if (length(rgb) != 3 || any(is.na(rgb)))
            stop("In Color mode, 'col' must be a valid color")
    } else {
        rgb = as.numeric(c(col, 0, 0))
        if (length(rgb) != 3 || any(is.na(rgb)))
            stop("In Grayscale mode, 'col' must be a scalar value")
    }

    n_row=dim(img)[2]
    n_col=dim(img)[1]
    boolv=vector(mode = "logical",n_row*n_col)

    boolv[x+n_col*y]<-TRUE
    if(size>0){
        for(i in 1:(size)){
            boolv[(x+i)+n_col*(y+i)]<-TRUE
            boolv[(x+i)+n_col*(y-i)]<-TRUE
            boolv[(x-i)+n_col*(y+i)]<-TRUE
            boolv[(x-i)+n_col*(y-i)]<-TRUE
        }
    }

    boolm=matrix(boolv,nrow=n_row,ncol=n_col)
    img[boolm]<-col
    invisible(img)
}

drawLine<-function (img, x1, y1, x2, y2, col=0.75, z = 1){
    #EBImage:::validImage(img)
    if (EBImage::colorMode(img) == EBImage::Color)
        stop("this method doesn't support the 'Color' color mode")
    if (any(is.na(img)))
        stop("'img' shouldn't contain any NAs")
    if (missing(x1)) stop("'x1' is missing")
    if (missing(y1)) stop("'y1' is missing")
    if (missing(x2)) stop("'x2' is missing")
    if (missing(y2)) stop("'y2' is missing")
    if (z < 1 | z > EBImage::numberOfFrames(img, "render"))
        stop("'z' must be a positive integer lower than the number of image frames")
    xy2z = as.integer(c(x1, y1, x2, y2, z - 1))
    if (length(xy2z) != 5 || any(is.na(xy2z)))
        stop("'x1', 'y1', 'x2', 'y2' and 'z' must be scalar values")
    if (EBImage::colorMode(img) == EBImage::Color) {
        rgb = as.numeric(col2rgb(col)/255)
        if (length(rgb) != 3 || any(is.na(rgb)))
            stop("In Color mode, 'col' must be a valid color")
    } else {
        rgb = as.numeric(c(col, 0, 0))
        if (length(rgb) != 3 || any(is.na(rgb)))
            stop("In Grayscale mode, 'col' must be a scalar value")
    }

    boolm=matrix(FALSE,nrow=dim(img)[1],ncol=dim(img)[2])
    if(x1==x2){
        boolm[x1,min(y1,y2):max(y1,y2)]<-TRUE
    }else if(y1==y2){
        boolm[min(x1,x2):max(x1,x2),y1]<-TRUE
    }else{
        for(i in min(x1,x2):max(x1,x2))
            for(j in min(y1,y2):max(y1,y2))
                boolm[i,j]=abs((y1-y2)*i+(x2-x1)*j+x1*y2-x2*y1)<abs(x1-x2)/sqrt(2)
    }
    img[boolm]<-col
    invisible(img)
}

#work arround to replace EBImage deprcated function drawtext
drawText <- function(img,
                     labels,
                     x = NULL,
                     y = NULL,
                     adj = c(0,0),
                     reuseLabels = TRUE,
                     col = NULL) {

    img.width <- dim(img)[1]
    img.height <- dim(img)[2]

    if(is.null(x)) x <- img.width / 2
    if(is.null(y)) y <- img.width / 2

    #create name to uniquely identify label temporary tif file
    img.fname <- paste(tempdir(),
                       "/Rcell_drawText_",
                       paste(digest::digest(img),
                             digest::digest(labels),
                             sum(x),
                             sum(y),
                             adj[1],
                             adj[2],
                             sep = "_"),
                       ".tif",sep="")

    #create file if it doesn't exist
    if (!file.exists(img.fname) | !reuseLabels){

        open.dev <- dev.list()

        #create device
        tiff(filename = img.fname,
             width = img.width,
             height = img.height,
             units = "px",
             pointsize = 12,
             bg = "white",
             family = "")

        # change graphical parameters
        op <- par(mar = c(0,0,0,0), oma = c(0,0,0,0))
        if(!is.null(open.dev)) on.exit(par(op))

        #plot current img in device
        EBImage::display(img,method = "raster")

        #add text to img
        text(x / img.width,
             1 - y / img.height,
             labels,
             adj = adj,
             col = col)

        #close device
        dev.off(base::setdiff(dev.list(), open.dev))
    }

    #open image and change colorMode
    ##capture.output to avoid anoying text
    capture.output(img <- EBImage::getFrame(EBImage::readImage(img.fname), 1))
    EBImage::colorMode(img) <- EBImage::Grayscale

    if (img.width != dim(img)[1] | img.height != dim(img)[2]) {
        stop("drawText changed img size")
        }

    return(img)
}

#####################Private Functions#####################################

#pform = pformula
#side an integer specifying which side of the plot the axis is to be drawn on.
#The axis is placed as follows: 1=below, 2=left, 3=above and 4=right.
#To plot de axis, a tiff device is opened in a temp dir, text and lines are called, and the image is loaded to a EBImage Image
if(getRversion() >= "2.15.1")
    utils::globalVariables(c("img_coord","pos","cellID"))

.axis <- function(imgdf,
                  pform,
                  side = 1,
                  img.size = 41,
                  border = 1,
                  font.size = 14,
                  max.nchar = floor(img.size * 2 / font.size),
                  mex = 1.2,
                  bg.col = "white"){ #,font.col="black",line.col=0

    #create name to uniquely identify label temporary tif file
    img.fname <- paste("Rcell_axis",
                       digest::digest(imgdf),
                       digest::digest(pform),
                       side,
                       img.size,
                       border,
                       font.size,
                       max.nchar,
                       mex,
                       sep = "_")

    img.fname <- paste(tempdir(), "/", img.fname, ".tif", sep = "")

    if(!file.exists(img.fname)) { #if axis not previously created, create them

        if (side == 1) { #bottom axis
            term = pform$rterm
            imgdf$img_coord <- imgdf$img_x
        } else if (side == 2) { #left
            term = pform$lterm
            imgdf$img_coord <- imgdf$img_y
        } else if (side == 3) { #top
            term = pform$rterm
            imgdf$img_coord <- imgdf$img_x
        } else if (side == 4) { #rigth
            term = pform$lterm
            imgdf$img_coord <- imgdf$img_y
        }

        max_coord <- max(imgdf[,"img_coord"])
        axis.img.label.x <- c()
        axis.img.label.y <- c()
        axis.img.label.text <- c()
        axis.img.line.x <- c()
        axis.img.line.y <- c()

        if("." %in% term | length(term) == 0) { #no axis
            axis.img.dim <- c((img.size + border) * max_coord + border, 1)
            # axis.img<-Image(bg.col,colormode="Grayscale", dim = axis.img.dim) #draw on EBImage
        } else { #plot axis
            if(!all(c(term,"img_coord") %in% names(imgdf))) stop("some variable in formula not in imgdf")
            axisdf <- plyr::ddply(subset(imgdf, select = c("img_coord", term)), plyr::.(img_coord), unique)
            n_term <- length(term)
            axis.img.dim <- c((img.size + border) * max_coord + border, round(n_term * font.size * mex))
            # axis.img <- Image(bg.col,colormode="Grayscale", dim = axis.img.dim) #draw on EBImage
            axisdf$max_coord <- axisdf[,"img_coord"]
            axisdf$min_coord <- axisdf[,"img_coord"]

            for(i in 1:n_term){
                if(i > 1) {
                    axisdf <- plyr::ddply(axisdf,term[1:(n_term - i + 1)],
                                    function(df) {
                                        data.frame(img_coord = mean(df$img_coord),
                                                   min_coord = min(df$min_coord),
                                                   max_coord = max(df$max_coord),
                                                   df[1, term[i:n_term]])
                                    })


                    # do.call(c,dlply(axisdf,.(img_coord),function(df)
                    # seq((df[1,"min_coord"]-1)*(img.size+border)+3*border,df[1,"max_coord"]*(img.size+border)-3*border)
                    # ))
                    # if(side %in% c(2,4)) xline<-rev(dim(axis.img)[1]-xline) #draw on EBImage

                    axis.img.line.x <- c(axis.img.line.x, #draw on device
                                         do.call(c,
                                                 plyr::dlply(axisdf,
                                                       plyr::.(img_coord),
                                                       function(df) {
                                                           c((df[1, "min_coord"] - 1) * (img.size + border) + 3 * border,
                                                             df[1, "max_coord"] * (img.size + border) - 3 * border, NA)
                                                       })
                                         )
                    )

                    ylinepos <- switch(side,
                                       round(mex * font.size * (i - 1)),	 	#side==1
                                       round((n_term - i + 1) * mex * font.size),	#side==2
                                       round((n_term - i + 1) * mex * font.size),	#side==3
                                       round(mex * font.size * (i - 1)))	#side==4

                    axis.img.line.y <- c(axis.img.line.y, #draw on device
                                         do.call(c, plyr::dlply(axisdf, plyr::.(img_coord),
                                                          function(df) {
                                                              c(ylinepos, ylinepos, NA)})))

                    #axis.img[as.integer(xline),ylinepos]<-line.col #draw on EBImage
                }

                lab = as.character(axisdf[,term[n_term - i + 1]])
                if(term[n_term-i+1]=="sample"){ #dealing with "sample" labels
                    if(all(c("pos","cellID","sample")%in%names(imgdf))){ #pos and cellID columns in imgdf
                        pIDs<-plyr::ddply(subset(imgdf,select=c("pos","cellID","sample")), plyr::.(pos,cellID)
                                    ,function(df)data.frame(sample=ifelse(length(unique(df$sample))==1,unique(df$sample),NA)))
                        pIDs<-transform(pIDs,lab=as.character(interaction(pIDs$pos,pIDs$cellID)))
                        dtsl<-dim(table(pIDs$sample,pIDs$lab))
                        if(any(is.na(pIDs[,3]))|(dtsl[2]>dtsl[1])){ #not unique mapping between sample number and pos.cellID
                            lab=rep(".",times=length(lab)) #using dots
                        }else { #unique mapping #using pos.cellID
                            lab=plyr::join(axisdf,pIDs,by="sample")$lab
                        }
                    } else { #NO pos and cellID columns
                        lab=rep(".",times=length(lab)) #using dots
                    }
                }

                lab=substr(lab,1,(axisdf$max_coord-axisdf$min_coord+1)*max.nchar)  #recorto labels para que entren
                #lab=substr(lab,1,max.nchar) #allow more flexbility here

                ylabelpos<-	switch(side
                                   ,round((1+(i-1)*mex)*font.size)		#side==1
                                   ,round((n_term-i+2-mex)*mex*font.size)	#side==2
                                   ,round((n_term-i+2-mex)*mex*font.size)	#side==3
                                   ,round((1+(i-1)*mex)*font.size)	)	#side==4
                xlabelpos<-	switch(side
                                   ,(axisdf$img_coord-0.5)*(img.size+border) - .nchar(lab)*floor(font.size/4) 				#side==1
                                   ,(max_coord-axisdf$img_coord+0.5)*(img.size+border) - .nchar(lab)*floor(font.size/4)		#side==2
                                   ,(axisdf$img_coord-0.5)*(img.size+border) - .nchar(lab)*floor(font.size/4) 				#side==3
                                   ,(max_coord-axisdf$img_coord+0.5)*(img.size+border) - .nchar(lab)*floor(font.size/4)	)	#side==4

                axis.img.label.x<-c(axis.img.label.x,xlabelpos) #draw on device
                axis.img.label.y<-c(axis.img.label.y,rep(ylabelpos,times=length(xlabelpos))) #draw on device
                axis.img.label.text<-c(axis.img.label.text,lab) #draw on device

                # axis.img<-drawText(axis.img,labels=lab,x=xlabelpos,y=ylabelpos) #draw on EBImage

            }
        }

        #drawing on device
        #create device
        tiff(filename = img.fname,
             width = axis.img.dim[1], height = axis.img.dim[2], units = "px", pointsize = 12,
             bg = "white", family = "")

        #change graphical parameters
        op<-par(mar=c(0,0,0,0),oma=c(0,0,0,0),xaxs="i",yaxs="i")
        on.exit(par(op))

        #new plot
        plot.new()

        #add text to device
        text((axis.img.label.x-1)/(axis.img.dim[1]-1),1-(axis.img.label.y-1)/(axis.img.dim[2]-1),axis.img.label.text,adj=c(0,0))

        #add lines to device
        if(side==1){
            lines((axis.img.line.x-1)/(axis.img.dim[1]-1),1-(axis.img.line.y-1)/(axis.img.dim[2]-1))
        }else if(side==2){
            lines(1-(axis.img.line.x-1)/(axis.img.dim[1]-1),1-(axis.img.line.y-1)/(axis.img.dim[2]-1))
        }

        #close device
        dev.off()
    }

    #open image and change colorMode
    capture.output(axis.img<-EBImage::getFrame(EBImage::readImage(img.fname),1))#capture.output to avoid anoying text
    EBImage::colorMode(axis.img)<-EBImage::Grayscale


    if(side %in% c(2,4)) axis.img<-EBImage::rotate(axis.img,filter="none",angle=(-90),bg.col=bg.col)
    #,output.dim=rev(dim(axis.img)))
    #,output.origin=c(0,dim(axis.img)[1]-dim(axis.img)[2]-1))
    return(axis.img)

}

if(getRversion() >= "2.15.1")
    utils::globalVariables(c("facet_z","facet_id"))
.append.panel.idxyz<-function(imgdf,pform,var.prefix="panel_",nx=NULL,ny=NULL,rev.y=TRUE){
    if(pform$type=="none"){
        tmp<-plyr::summarise(imgdf,facet_id=factor(1),facet_x=1,facet_y=1)
        names(tmp)<-paste(var.prefix,c("id","x","y"),sep="")
        imgdf<-cbind(imgdf,tmp)
    } else if (pform$type=="wrap_horizontal"){
        imgdf<-cbind(imgdf,
                     .wrap(subset(imgdf,select=pform$rterm),nx=nx,ny=ny,vertical=FALSE,var.prefix=var.prefix))
    } else if (pform$type=="wrap_vertical"){
        imgdf<-cbind(imgdf,
                     .wrap(subset(imgdf,select=pform$lterm),nx=nx,ny=ny,vertical=TRUE,var.prefix=var.prefix))
    } else if (pform$type=="grid"){
        imgdf<-cbind(imgdf,
                     .grid(imgdf,pform$lterm,pform$rterm,var.prefix=var.prefix,rev.y=rev.y))
    } else stop("unknown formula type")

    if(var.prefix=="img_") return(imgdf)

    #dealing with slice (z) term
    if(length(pform$sterm)>0){ #multi layer image
        imgdf$facet_z <- as.numeric(.ordered.interaction(subset(imgdf,select=pform$sterm)))
        imgdf$facet_id <- .ordered.interaction(subset(imgdf,select=c(facet_z,facet_id)))
    } else { #single layer image
        imgdf<-transform(imgdf,facet_z=1)
    }

    return(imgdf)
}

if(getRversion() >= "2.15.1")
    utils::globalVariables(c("factor_id_x","factor_id_y","factor_id","factor_x","factor_y"))

.grid <- function(df,lterm,rterm,var.prefix="panel_",rev.y=TRUE){

    df$factor_id_x<-.ordered.interaction(subset(df,select=rterm))
    df$factor_id_y<-.ordered.interaction(subset(df,select=lterm))
    df$factor_id<-.ordered.interaction(subset(df,select=c(factor_id_x,factor_id_y)),sep="_._")
    df$factor_x<-as.numeric(df$factor_id_x)
    df$factor_y<-as.numeric(df$factor_id_y)

    if(isTRUE(var.prefix=="img_")&isTRUE(rev.y))
        df$factor_y<-max(df$factor_y)-df$factor_y+1 #inverting y axis for imgs

    df<-subset(df,select=c(factor_id,factor_id_x,factor_id_y,factor_x,factor_y))

    names(df)<-paste(var.prefix,c("id","id_x","id_y","x","y"),sep="")

    return(df)
}

if(getRversion() >= "2.15.1") utils::globalVariables(c("factor_id"))
.wrap <- function(df,nx=NULL,ny=NULL,vertical=FALSE,var.prefix="panel_"){

    df$factor_id <- .ordered.interaction(df)

    ntot=nlevels(df$factor_id)
    if(is.null(nx)&is.null(ny)){
        nx=ceiling(sqrt(ntot))
        ny=ceiling(ntot/nx)
    } else if(is.null(nx)){
        nx=ceiling(ntot/ny)
    } else if(is.null(ny)){
        ny=ceiling(ntot/nx)
    } else if(nx*ny<ntot) stop("Not enough space for all panels. Increase nx or ny.")

    if(!vertical){ #horizontal wrapping
        idxy<-data.frame(levels(df$factor_id),(0:(ntot-1)%%nx)+1,(0:(ntot-1)%/%nx)+1)
    } else { #horizontal vertical
        idxy<-data.frame(levels(df$factor_id),(0:(ntot-1)%/%ny)+1,(0:(ntot-1)%%ny)+1)
    }
    names(idxy)<-c("factor_id","factor_x","factor_y")
    df<-plyr::join(subset(df,select=factor_id),idxy,by="factor_id")
    names(df)<-paste(var.prefix,c("id","x","y"),sep="")
    return(df)
}

.ordered.interaction <- function(df,drop=TRUE,sep="_"){
    df.names<-names(df)
    df$ord.int<-interaction(df,drop=drop,sep=sep)
    return(
        factor(df$ord.int,ordered=TRUE
               ,levels=reshape::sort_df(plyr::ddply(df,df.names,function(df)df[1,]),df.names)$ord.int)
    )
}

if(getRversion() >= "2.15.1") utils::globalVariables(c("ucid"))
.sample.N.ucid<-function(df,N){

    group.ucid=unique(df$ucid)
    len.group=length(group.ucid)
    if(is.null(N)) N=len.group
    if(len.group==0) stop("Empty group selected",call. = FALSE)
    if(len.group==1){
        df<-transform(df,sample=1)
    }else if(len.group<=N){ #less cells in group than required
        df<-plyr::join(df,data.frame(ucid=sample(group.ucid,length(group.ucid)),sample=1:len.group),by="ucid")
    } else { #more cells than required
        group.ucid<-sample(group.ucid,N)
        df<-subset(df,ucid %in% group.ucid)
        df<-plyr::join(df,data.frame(ucid=group.ucid,sample=1:N),by="ucid")
    }
    return(df)
}

.parse.formula <- function(formula,force.valid.names=FALSE){

    if (length(formula) == 3) { # right and left term

        # Rcell::.parseSide
        lterm <- .parseSide(formula[[2]])

    } else if(length(formula) == 2) { #only right term
        lterm <- "."
    }

    rterm <- formula[[length(formula)]]
    if (length(rterm) == 3 && rterm[[1]] == as.name("|")){ #slice term
        sterm <- .parseSide(rterm[[3]])
        rterm <- .parseSide(rterm[[2]])
    } else { #no slice term
        sterm <- c()
        rterm <- .parseSide(rterm)
    }

    #checking that "." is not used with other vars
    if(length(lterm) > 1 & is.element(".", lterm)) {
        stop("left term of formula has . and other variables")}

    if(length(rterm) > 1 & is.element(".", rterm)) {
        stop("right term of formula has . and other variables")}

    if(is.element(".", sterm)){
        stop("dot not allowed in slice term")}

    #dealing with "sample" synomins
    lterm[lterm %in% c("...", "cell", "cells")] <- "sample"
    rterm[rterm %in% c("...", "cell", "cells")] <- "sample"
    sterm[sterm %in% c("...", "cell", "cells")] <- "sample"

    if("." %in% lterm & "." %in% rterm) {
        type = "none"
    } else if("." %in% lterm) {
        type = "wrap_horizontal"
    } else if("." %in% rterm) {
        type = "wrap_vertical"
    } else {
        type="grid"
    }

    pform <- list(lterm = lterm, rterm = rterm, sterm = sterm, type = type)

    if(isTRUE(force.valid.names))
        pform <- lapply(pform, make.names)

    return(pform)
}


#auxiliary functions for .parse.formula
.parseSide <- function(model) {
    model.vars <- list()
    while (length(model) == 3 && model[[1]] == as.name("+")) {
        model.vars <- c(model.vars, model[[3]])
        model <- model[[2]]
    }
    return(sapply(rev(c(model.vars, model)), .expr2char))
}


.expr2char <- function(x) paste(deparse(x), collapse = "")

###code adapted from check_formula (reshape package) from Hadley Wickham
# Check formula
# Checks that formula is a valid reshaping formula.
#
# \enumerate{
#   \item variable names not found in molten data
#   \item same variable used in multiple places
# }
# @arguments formula to check
# @arguments vector of variable names
# @keyword internal

.check_formula <- function(formula, varnames) {
    vars <- all.vars(formula)
    unknown <- base::setdiff(vars, c(".", "...", "cell", "cells", varnames))

    if (length(unknown) > 0) {
        stop("formula contains variables not found in cell.image img.desc: ",
             paste(unknown, collapse = ", "),
             call. = FALSE)
    }

    vars <- vars[vars != "."]

    if (length(unique(vars)) < length(vars)) {
        stop("Variable names repeated", call. = FALSE)
    }
}

.unique.na <- function(x){
    ux <- unique(x)
    if(length(ux) == 1) {
        return(ux)
    } else {
        return(NA)
    }
}
