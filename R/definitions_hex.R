
#' Funcion para armar un dataframe para hacer plots tipo HEX
#'
#' @param cdata A "cdata" Rcell data.frame (not the object).
#' @param facetVars are "pdata" variables
#' @param varx the x-axis "pdata" variable
#' @param vary the y-axis "pdata" variable
#' @return Lots of stuff.
#' @import hexbin
hexPlotDf<- function(cdata, facetVars = c("pos", "treatment"), varx = "a.tot", vary = "el.p"){
  print("F9")
  # Reemplacé las funciones anteriores por esta, que procesa un cdata y lo pone lindo para geom_hex
  # Lindo según: https://stackoverflow.com/questions/14495111/setting-hex-bins-in-ggplot2-to-same-size

  if(length(facetVars) != 0){
    # Hacer un dataframe con un factor que combine todos los facets
    dfactor <- cdata[,c(facetVars), drop = FALSE]
    cdata$factor <- as.factor(
      apply(dfactor, 1, function(x) paste(x, collapse = "_"))
    )
    dfactor <- unique(cdata[,c("factor", facetVars)])

    H <- myhexbin(cdata, varx, vary, xbins = 25)

    # Siguiendo: https://stackanswers.net/questions/setting-hex-bins-in-ggplot2-to-same-size
    counts <- hexbin::hexTapply(H, cdata$factor, table)
    counts <- t(simplify2array(counts))
    counts <- reshape2::melt(counts)
    colnames(counts) <- c("cID", "factor", "counts")

    h <- data.frame(hcell2xy(H), cID = H@cell)
    h <- merge(counts, h)
    h$counts[h$counts == 0] <- NA  # Para que después los valores NA no se vean (serán transparentes en ggplot)

    h <- merge(h, dfactor, by = "factor")
    h <- dplyr::select(h, -factor)
  } else {
    # No considerar facets
    # cdata2 <- cdata
    # facetVars = c()
    H <- myhexbin(cdata, varx, vary, xbins = 25)
    h <- data.frame(hcell2xy(H), cID = H@cell, counts = H@count)
  }

  #axisRatio <- (max(h$x) - min(h$x))/(max(h$y) - min(h$y))

  return(h)

  # ggplot(h, aes(x=x, y=y, fill = counts)) +
  #     geom_hex(stat="identity") +
  #     coord_equal (axisRatio) +
  #     theme_bw()+ xlab("a.tot") + ylab("el.p") +
  #     scale_fill_continuous (low = "grey80", high = "red", na.value = "#00000000") +
  #     facet_grid(pH~treatment)
}

#' Función hexbin
#' Esencialmente hexbin() pero un poco más a mi gusto
myhexbin <- function(bindata, varx , vary, xbins = 25){
  print("F10")

  # Ranges are single valued if drawing only one polygon, fixed here:
  bindata <- as.data.frame(bindata)
  if(nrow(bindata) == 1) {
      xbnds <- c(bindata[,varx]*0.99, bindata[,varx]*1.01)
    } else {
      xbnds <- bindata[,varx]
    }

  if(nrow(bindata) == 1) {
      ybnds <- c(bindata[,vary]*0.99, bindata[,vary]*1.01)
    } else {
      ybnds <- bindata[,vary]
    }

  h <- hexbin::hexbin(bindata[,varx], bindata[,vary], xbins = xbins, IDs = TRUE,
                      xbnds = range (xbnds),
                      ybnds = range (ybnds))
  return(h)
}
