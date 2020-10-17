#' Hu moments
#'
#' @param xy Matriz de dos columnas (dim1 y dim2), con las coordenadas XY del boundary binario. No es importante el orden de los puntos.
#' @return Los 7 Hu moments + 1 extra (?)
# @example
# xy <- EBImage::imageData(readImage("inst/hu_moments_img2.png"))
# xy <- which(xy>0, arr.ind = TRUE)
# plot(xy)
# xy.hu <- hu.moments(xy)
# xy.hu <- -sign(xy.hu)*log10(abs(xy.hu))
hu.moments2 <- function(xy){
  ## Define function for calculating eta
  eta <- function(xy, i, j){
    ## Calculate central moment mu(i,j)
    mu <- sum(((x-cent[1])**i)*((y-cent[2])**j)*px.intensity)
    
    ## Calculate and return normalized central moment eta(i,j)
    mu/(mu00**((i+j)/2+1))
  }
  
  px.intensity <- 255
  
  x <- xy[,1]
  y <- xy[,2]
  
  ## Calculate basic object moments M(0,0) M(0,1) M(1,0) M(1,1)
  M <- matrix(sapply(1:4, function(n,i,j,x,y) sum((x**i[n])*(y**j[n])*px.intensity),
                     x=x, y=y, i=c(0,1,0,1), j=c(0,0,1,1)),2)
  
  ## Calculate object centroid (x,y)
  cent <- diag(apply(M/M[1,1],2,rev))
  
  ## Calculate object mu(0,0)
  mu00 <- sum(((x-cent[1])**0)*((y-cent[2])**0)*px.intensity)
  
  I1 <- eta(xy,2,0) + eta(xy,0,2)
  I2 <- (eta(xy,2,0) - eta(xy,0,2))**2 + 4*eta(xy,1,1)**2
  I3 <- (eta(xy,3,0) - 3*eta(xy,1,2))**2 + (3*eta(xy,2,1) - eta(xy,0,3))**2 #not very useful?
  I4 <- (eta(xy,3,0) + eta(xy,1,2))**2 + (eta(xy,2,1) + eta(xy,0,3))**2
  I5 <- (eta(xy,3,0) - 3*eta(xy,1,2))*(eta(xy,3,0) + eta(xy,1,2))*
    ((eta(xy,3,0) + eta(xy,1,2))**2 - 3*(eta(xy,2,1) + eta(xy,0,3))**2) +
    (3*eta(xy,2,1) - eta(xy,0,3)) * (eta(xy,2,1) + eta(xy,0,3)) *
    (3*(eta(xy,3,0) + eta(xy,1,2))**2 - (eta(xy,2,1) + eta(xy,0,3))**2)
  I6 <- (eta(xy,2,0) - eta(xy,0,2)) * ((eta(xy,3,0) + eta(xy,1,2))**2 -
                                         (eta(xy,2,1) + eta(xy,0,3))**2) +
    4*eta(xy,1,1) * (eta(xy,3,0) + eta(xy,1,2)) * (eta(xy,2,1) + eta(xy,0,3))
  I7 <- (3*eta(xy,2,1) - eta(xy,0,3)) * (eta(xy,3,0) + eta(xy,1,2)) *
    ((eta(xy,3,0) + eta(xy,1,2))**2 - 3*(eta(xy,2,1) + eta(xy,0,3))**2) -
    (eta(xy,3,0) - 3*eta(xy,1,2)) * (eta(xy,2,1) + eta(xy,0,3)) *
    (3*(eta(xy,3,0) + eta(xy,1,2))**2 - (eta(xy,2,1) + eta(xy,0,3))**2)
  I8 <- eta(xy,1,1) * ((eta(xy,3,0) + eta(xy,1,2))**2 - (eta(xy,0,3) + eta(xy,2,1))**2) -
    (eta(xy,2,0) - eta(xy,0,2)) * (eta(xy,3,0) + eta(xy,1,2)) * (eta(xy,0,3) + eta(xy,2,1))
  
  #hu <- c(I1, I2, I3, I4, I5, I6, I7, I8)
  #return(-sign(hu)*log10(abs(hu)))
  return(c(I1, I2, I3, I4, I5, I6, I7, I8))
}
