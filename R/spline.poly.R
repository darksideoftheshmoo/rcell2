#' whuber's spline.poly as shared in StackExchange
#' 
#' See: https://gis.stackexchange.com/a/24929
#'
#' Splining a polygon with \code{stats::spline}
#'
#'   The rows of 'xy' give coordinates of the boundary vertices, in order.
#'   'vertices' is the number of spline vertices to create.
#'              (Not all are used: some are clipped from the ends.)
#'   'k' is the number of points to wrap around the ends to obtain
#'       a smooth periodic spline.
#'
#' @return Returns an array of points. 
#' @export
#' @param xy A data.frame with ordered X/Y pairs representing the polygon
#' @param vertices The amount of vertices in the final smoothing/interpolation.
#' @param k The amount of points used to "close" the path, and prevent discontinuous derivatives at the end-points.
#' 
spline.poly <- function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- stats::spline(x=1:(n+2*k),
                               y=data[,1], 
                               n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}

#' whuber's spline.poly adapted to stats::smooth.spline
#' 
#' See: https://gis.stackexchange.com/a/24929
#'
#' Splining a polygon with \code{stats::smooth.spline}
#'
#'   The rows of 'xy' give coordinates of the boundary vertices, in order.
#'   'vertices' is the number of spline vertices to create.
#'              (Not all are used: some are clipped from the ends.)
#'   'k' is the number of points to wrap around the ends to obtain
#'       a smooth periodic spline.
#'
#' @return Returns a data.frame of points. 
#' @export
#' @param xy A data.frame with ordered X/Y pairs representing the polygon
#' @param k The amount of points used to "close" the path, and prevent discontinuous derivatives at the end-points.
#' @param dof The degrees of freedom used by \code{stats::smooth.spline}
#'
smooth.spline.poly <- function(xy, k=3, dof=5, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- stats::smooth.spline(x=1:(n+2*k),
                                      y=data[,1,drop=T], 
                                      df=dof)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- stats::smooth.spline(
    x=1:(n+2*k), 
    y=data[,2,drop=T],
    df=dof
  )$y
  
  # Retain only the middle part.
  result <- cbind(x = x1, y = x2)[k < x & x <= n+k, ]
  
  result <- as.data.frame(result)
  
  result
}
