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
  n_vertices <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n_vertices-k+1):n_vertices,], # add k-last points to the beginning
                  xy, 
                  xy[1:k, ] # add k-first points to the end
                  )
  } else {
    data <- xy
  }
  
  # indices vector (for k-wraped data)
  x.idxs <- 1:(n_vertices+2*k)
  
  # Spline the x and y coordinates:
  data.spline.1 <- stats::smooth.spline(x=x.idxs,
                                        y=data[,1,drop=T], 
                                        df=dof)
  
  data.spline.2 <- stats::smooth.spline(x=x.idxs, 
                                        y=data[,2,drop=T],
                                        df=dof)
  
  # NEW: interpolate output to constant length ###
  # Prepare re-sampled/re-interpolated results
  # for constant length output == length(x.idxs.new)
  # x.idxs.new <- seq(from=1, to=max(x.idxs), 
  x.idxs.new <- seq(from=k, # this excludes indices < k
                    # using "k" instead of "k+1" closes the path, leaving no gap
                    to=length(x.idxs)-k, # excludes k-points added to the end
                    length.out=length.out+1)
  x1 <- predict(data.spline.1, x=x.idxs.new)$y
  x2 <- predict(data.spline.2, x=x.idxs.new)$y
  
  # Get x value ranges for (unwraped) row indices
  unwrap.filter <- (x.idxs.new > k) & (x.idxs.new <= n_vertices+k)
  
  # Prepare results, retaining only the middle part
  result <- cbind(x = x1, y = x2)[unwrap.filter,]
  
  # Convert to data.frame
  result <- as.data.frame(result)
  
  return(result)
}
