#' Splining a polygon with \code{stats::spline}
#' 
#' whuber's spline.poly as shared in StackExchange
#' 
#' See: https://gis.stackexchange.com/a/24929
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
#' @return Returns a data.frame of smooth points. 
#' @export
#' @param xy A data.frame with ordered X/Y pairs representing the polygon
#' @param k The amount of points used to "close" the path, and prevent discontinuous derivatives at the end-points.
#' @param dof The degrees of freedom used by \code{stats::smooth.spline}
#' @param length.out Set to an integer length for the smoothed output (defaults to input length).
#'
smooth.spline.poly <- function(xy, k=3, dof=5, length.out=nrow(xy), ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # indices vector
  x.idxs <- 1:(n+2*k)
  
  # Spline the x and y coordinates:
  data.spline.1 <- stats::smooth.spline(x=x.idxs,
                                        y=data[,1,drop=T], 
                                        df=dof)
  
  data.spline.2 <- stats::smooth.spline(x=x.idxs, 
                                        y=data[,2,drop=T],
                                        df=dof)
  
  # NEW: interpolate output to constant length ###
  
  # Get x value ranges for (unwraped) row indices
  x.range <- range(x.idxs[unwrap.filter])
  # Interpolate them to a new length
  x.idxs.new <- seq(from=x.range[1], to=x.range[2], length.out=length.out)
  
  # Prepare re-sampled/re-interpolated results
  # for constant length output == length(x.idxs.new)
  x <- x.idxs.new
  x1 <- predict(data.spline.1, x = x.idxs.new)$y
  x2 <- predict(data.spline.2, x = x.idxs.new)$y
  
  # Prepare results
  result <- cbind(x = x1, y = x2)
  
  # OLD: output has variable length; i.e. it depends on input length ###
  
  # Prepare results
  # x <- data.spline.1$x
  # x1 <- data.spline.1$y
  # x2 <- data.spline.2$y
  
  # Retain only the middle part
  # unwrap.filter <- (k < x) & (x <= n+k)
  # result <- cbind(x = x1, y = x2)[unwrap.filter, ]
  
  # Convert to data.frame
  result <- as.data.frame(result)
  
  return(result)
}
