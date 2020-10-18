## Moments
xy <- mask.data[[1]]$cell.px
xy <- imageData(readImage('/Users/andreas/Downloads/K0_2.tif'))
xy2 <- xy
xy <- which(xy>0, arr.ind = TRUE)

## Arbitrary moment M(i,j)
mment <- function(xy, i, j){
  x <- xy[,1]
  y <- xy[,2]
  sum((x**i)*(y**j)*255)
}
mment(xy, 0, 0)
mment(xy, 0, 1)
mment(xy, 1, 0)
mment(xy, 1, 1)

## Arbitrary moment 2 M(i,j)
mment2 <- function(xy, i, j){
  x <- 1:dim(xy)[1]
  y <- 1:dim(xy)[2]
  xyprod <- (x**i)%*%t(y**j)
  print(sum(xyprod * xy * 255))
}
mment2(xy2, 0, 0)
mment2(xy2, 1, 0)
mment2(xy2, 0, 1)
mment2(xy2, 1, 1)

## Moments (0,0) (1,0) (0,1) (1,1)
tx <- proc.time()
M <- matrix(sapply(1:4, function(xy, i, j, n) mment(xy,i[n],j[n]), xy=xy, i=c(0,1,0,1), j=c(0,0,1,1)),2)
print(proc.time()-tx)

## Centroid 
M.centroid <- diag(apply(M/M[1,1],2,rev))
#xbar <- M[2,1]/M[1,1]
#ybar <- M[1,2]/M[1,1]

## Central moments (translation invariant)
MU <- matrix(sapply(1:4, function(n,i,j,x,y,cent){
  sum(((x-cent[1])**i[n])*((y-cent[2])**j[n])*255)
},x=xy[,1], y=xy[,2], i=c(0,1,0,1), j=c(0,0,1,1), cent=M.centroid),2)

## Arbitrary central moment (i,j)
mu <- function(xy, i, j){
  x <- xy[,1]
  y <- xy[,2]
  
  M <- matrix(sapply(1:4, function(xy,i,j,n) mment(xy,i[n],j[n])
                     , xy=xy, i=c(0,1,0,1), j=c(0,0,1,1)),2)
  
  #M <- matrix(sapply(1:4, function(n,i,j,x,y) sum((x**i[n])%*%t(y)**j[n]),
  #                   x=x, y=y, i=c(0,1,0,1), j=c(0,0,1,1)),2)
  M.cent <- diag(apply(M/M[1,1],2,rev))
  
  sum(((x-M.cent[1])**i)*((y-M.cent[2])**j)*255)
  #sum(((x-M.cent[1])**i)%*%t(y-M.cent[2])**j)
}
mu(xy, 1, 1)

## Central centroid
MU.centroid <- diag(apply(MU/MU[1,1],2,rev))

## Normalized central moments (translation and scale invariant)
ETA <- matrix(sapply(1:4, function(n,i,j,M){
  M[i[n]+1,j[n]+1]/(M[1,1]**((i[n]+j[n])/2+1))
}, i=c(0,1,0,1), j=c(0,0,1,1), M=MU),2)


## Arbitrary normalized central moments
eta <- function(xy, i, j){
  x <- xy[,1]
  y <- xy[,2]
  
  intensity.px <- 255
  
  ## Calculate basic object moments: (0,0) (0,1) (1,0) (1,1)
  #M <- matrix(sapply(1:4, function(n,i,j,x,y) sum((x**i[n])*(y**j[n])*255),
  #                   x=x, y=y, i=c(0,1,0,1), j=c(0,0,1,1)),2)
  
  ## Calculate object centroid (x,y)
  #cent <- diag(apply(M/M[1,1],2,rev))
  
  ## Calculate central moment (i,j) and (0,0)
  mu <- sum(((x-cent[1])**i)*((y-cent[2])**j)*intensity.px)
  #mu00 <- sum(((x-cent[1])**0)*((y-cent[2])**0)*255)
  
  ## Calculate and return normalized central moment (j,j)
  mu/(mu00**((i+j)/2+1))
}
#eta(xy, 2, 2)

## Hu moments
hu.moments <- function(xy){
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

#xy <- mask.data[[24]]$cell.px
xy <- imageData(readImage('/Users/andreas/Downloads/K0_2.tif'))
xy <- imageData(readImage('/Users/andreas/Downloads/K0_3.tif'))
xy <- imageData(readImage('/Users/andreas/Downloads/K0_4.tif'))
xy <- imageData(readImage('/Users/andreas/Downloads/K0_5.tif'))
xy <- which(xy>0, arr.ind = TRUE)
plot(xy)

xy.hu <- hu.moments(xy)
xy.hu <- -sign(xy.hu)*log10(abs(xy.hu))
