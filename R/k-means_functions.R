#' X.normalize
#'
#' Normalize X
#'
#' @param X a data.frame with fit.vars
#'
#' @return Normalized X
#'
X.normalize <- function(X){
  X <- t((t(X)-apply(X, 2, mean))/apply(X, 2, sd))
  rownames(X) <- NULL
  colnames(X) <- NULL
  return(X)
}


#' add_fit_variables
#'
#' Create custom fit variables
#'
#' @param x a cell.data data.frame
#'
#' @return a cell.data data.frame to which new variables have been added
#'
add_fit_variables <- function(x){
  x <- transform(x, perim.bud=(pi*maj.axis))
  x <- transform(x, perim.hglass=(pi*min.axis))
  x <- transform(x, perim.circle=(pi*(min.axis/2)**2))
  
  x <- transform(x, a.bud=((pi/4)*(2*min.axis**2-2*min.axis*maj.axis+maj.axis**2)))
  x <- transform(x, a.hglass=(2*pi*(min.axis/2)**2))
  x <- transform(x, a.ellipse=(pi*min.axis*maj.axis))
  x <- transform(x, a.circle=(pi*(min.axis/2)**2))
  
  x <- transform(x, elongatedness=(maj.axis**2/min.axis))
  x <- transform(x, eccentricity=sqrt(1-((min.axis/2)/(maj.axis/2))**2))
  x <- transform(x, flattening=(maj.axis-min.axis)/(maj.axis+min.axis))
  x <- transform(x, perim.by.a=(perim/a.tot))
  x <- transform(x, perim.by.perimellipse=(perim/ellipse.perim))
  x <- transform(x, perim.by.perimbud=(perim/perim.bud))
  x <- transform(x, perim.by.spherevol=(perim/sphere.vol))
  
  x <- transform(x, a.by.abud=(a.tot/a.bud))
  x <- transform(x, surf.by.vol=a.surf/sphere.vol)
  x <- transform(x, surf.by.perim=a.surf/perim)
  x <- transform(x, fit.var1=fft.stat/elongatedness)
  x <- transform(x, fit.var2=fft.stat/a.surf)
  x <- transform(x, fit.var3=fft.stat/a.tot)
  x <- transform(x, fit.var4=fft.stat/sphere.vol)
  x <- transform(x, fit.var5=a.surf/elongatedness)
  x <- transform(x, fit.var6=sphere.vol/elongatedness)
  x
}


#' get.fit.variables
#'
#' Extract columns for clustering
#'
#' @param x 
#' @param vars "m" for morpho.vars, "f" for fluor.vars
#'
#' @return data.frame with fit.variables ("extr.vars") only
#'
get.fit.variables <- function(x, vars='m'){
  f.channels <- substr(unlist(regmatches(names(x), gregexpr(paste0("(f.tot.)+([a-z]{1}$)"), names(x)))), 7, 8)
  
  #id.vars <- c('pos', 'cellID', 'time', 'ucid', 't.frame')
  id.vars <- c('ucid','t.frame')
  
  #morpho.vars <- c('xpos', 'ypos', 'fft.stat', 'perim', 'maj.axis', 'min.axis', 'rot.vol', 'con.vol', 'a.surf', 'sphere.vol')
  morpho.vars <- c('a.tot','fft.stat', 'perim', 'maj.axis', 'min.axis', 'a.surf', 'sphere.vol', 'ellipse.perim')
  
  #area.vars <- c('a.tot', 'a.tot.p1', 'a.tot.m1', 'a.tot.m2', 'a.tot.m3')
  #area.vars <- c('a.tot')
  
  #f.vars <- c('f.tot', 'f.tot.p1', 'f.tot.m1', 'f.tot.m2', 'f.tot.m3')
  f.vars <- c('f.tot')
  f.vars <- sapply(f.vars, function(x, y) paste(x,y,sep='.'), y=f.channels)
  f.vars <- c(f.vars, 'f.c', 'f.r', 'f.y')
  
  #bg.vars <- c('a.local.bg', 'a.local', 'a.local2.bg', 'a.local2')
  bg.vars <- c()
  
  #bg.f.vars <- c('f.bg', 'f.local.bg', 'f.local2.bg')
  #bg.f.vars <- c('f.bg', 'f.local.bg')
  #bg.f.vars <- sapply(bg.f.vars, function(x, y) paste(x,y,sep='.'), y=f.channels)
  bg.f.vars <- c()
  
  default.vars <- names(x)
  
  x <- add_fit_variables(x)
  
  new.vars <- names(x)
  new.vars <- new.vars[!(new.vars%in%default.vars)]
  
  all.vars <- names(x)
  
  ## Defines whether fluorescence ('f') and/or morphology ('m') variables should be extracted
  extr.vars <- c()
  if(sum(grep("m",vars))){
    morpho.vars <- c(id.vars, morpho.vars, new.vars)
    extr.vars <- c(extr.vars,morpho.vars)
  }
  if(sum(grep("f",vars))){
    fluor.vars <- c(f.vars, bg.vars, bg.f.vars)
    extr.vars <- c(extr.vars,fluor.vars)
  }
  
  return(x[,names(x) %in% extr.vars])
}

#' kmeans_clustering
#'
#' Perform k-means clustering on Cell-ID data.
#'
#' @param x cell.data object
#' @param k non-negative integer. The desired number of clusters.
#' @param centroids optional integer vector of length 'k' defining the rows in 'x$data' to be chosen as starting centroids. If not provided, these are selected randomly.
#' @param max_iter non-negative integer. The maximal number of iterations allowed, if the algorithm has not yet converged.
#' @param resume logical. If TRUE the algorithm resumes clustering from previously assigned clusters.
#' @param vars optional character vector
#' 
#' @example 
#' 
#' # Define custom classification variables (optional; selected morphology variables used by default)
#' class.variables = NULL
#' #class.variables <- c('a.tot','ellipse.perim','perim','maj.axis','min.axis','sphere.vol')
#' 
#' # Set number k-mean classes
#' k_count <- 10
#' 
#' # Set max number of iterations (*)
#' max_iter <- 50 
#' 
#' # (*) The function kmeans_clustering runs the algorithm iteratively,
#' # and keeps track of how many cells changed class between each iteration.
#' # If for any iteration there were no class updates, the function considers
#' # the classification complete, and breaks the run. Otherwise, it continues
#' # the classification until "max.iter" iterations have been completed.
#' 
#' # Run k-means classification; returns Cell-ID object with column "k",
#' # indicating the assigned class for each row/cell
#' X <- kmeans_clustering(X, k=k_count, max_iter=max_iter, vars=class.variables)
#'
#' @return the original cell.data object with appended 'k' and 'k_dist' columns indicating the assigned class and Euclidean distance the assigned centroid, respectively.
#' @export
#'
kmeans_clustering <- function(x, k=10, max_iter=100, vars=NULL, resume=FALSE){
  if(!resume){
    ## Remove pre-existing 'k' and 'k_dist' columns
    x$data <- x$data %>% mutate(k = NULL, k_dist = NULL)
  }
  
  ## Get data to classify
  cdata <-x$data[x$data$qc,]
  #cellk <- subset(x$data, qc==TRUE)
  
  ## Append training variables
  cdata <- get.fit.variables(cdata)
  
  if(!is.null(vars)){
    vars <- c(c('ucid','t.frame','tag.type'),vars)
    cdata <- cdata[,names(cdata)%in%vars]
  }
  
  ## Remove rows with NAs
  n.cols <- ncol(cdata)
  na.idx <- unique(unlist(sapply(1:n.cols,function(x,i) which(is.na(x[,i])), x=cdata)))
  if(length(na.idx)) cdata <- cdata[-na.idx,]
  
  ## Extract ucids (id) and data (Xd)
  fit.vars <- names(cdata)[!(names(cdata)%in%c('ucid','t.frame','tag.type'))]
  id <- cdata[,c('ucid','t.frame')]
  Xd <- cdata[,fit.vars]
  
  ### RUN k-means classification
  ## Normalize X
  Xd <- X.normalize(Xd)
  
  ## Randomly pick k samples as initial centroids
  k.sample <- sample(dim(Xd)[1],k)
  k.means <- Xd[k.sample,]
  
  ## Set random initial centroid assignments
  k.labs.old <- sample(1:k,length(id),replace=TRUE)
  
  ## Loop over centroid calculation as long as labels change
  ## Limit loop to 100 iterations
  i.diff <- 1
  i.count <- 0
  
  while(i.diff>0 && i.count<max_iter){
    ## Calculate distances of each row in Xd to each centroid in k.means
    dists <- apply(k.means,1,function(x,k) sqrt(rowSums(sweep(x, 2, k)**2)),x=Xd)
    
    ## Find nearest centroid for each row in Xd
    k.labs <- apply(dists,1,function(x) which(x==min(x))[1])
    i.diff <- sum(k.labs!=k.labs.old)
    i.count <- i.count + 1
    message(paste0("Iteration #",i.count,". Updates: ",i.diff))
          
    ## Calculate new centroids
    k.means <- t(sapply(1:k, function(i,x,n) colMeans(x[which(n==i),]), x=Xd, n=k.labs))
    k.labs.old <- k.labs
  }
  
  ## Sum up data by ucid, centroid, and distance to centroid
  k.iter <- seq_len(k)
  dists <- apply(k.means,1,function(x,k) sqrt(rowSums(sweep(x, 2, k)**2)),x=Xd)
  k.idx <- lapply(k.iter,function(i,l) which(l==i), l=k.labs)
  k.dists <- lapply(k.iter,function(d,k,i) d[k[[i]],i],d=dists,k=k.idx)
  k.ucids <- lapply(k.iter,function(u,k,i) u[k[[i]],], u=id, k=k.idx)
  k.data <- mapply(function(ucid,k,k.dist) cbind(ucid,k,k.dist),k.dist=k.dists,ucid=k.ucids,k=k.iter,SIMPLIFY=FALSE)
  k.data <- lapply(k.data,function(x) x[order(x[,3]),])
  k.data <- do.call(rbind.data.frame, k.data)
  
  x$data <- merge(x$data, k.data, all.x=TRUE)
  
  return(x)
}