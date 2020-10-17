#### X.normalize ####
## Normalize X
X.normalize <- function(X){
  X <- t((t(X)-apply(X, 2, mean))/apply(X, 2, sd))
  rownames(X) <- NULL
  colnames(X) <- NULL
  return(X)
}


#### make.fit.variables ####
#Create custom fit variables; 
make.fit.variables <- function(x){
  x <- transform(x, perim.bud=(pi*maj.axis))
  x <- transform(x, perim.hglass=(pi*min.axis))
  x <- transform(x, perim.ellipse=(pi*(3*(maj.axis/2+min.axis/2)-sqrt((3*maj.axis/2+min.axis/2)*(maj.axis/2+3*min.axis/2)))))
  x <- transform(x, perim.circle=(pi*(min.axis/2)**2))
  
  x <- transform(x, a.bud=((pi/4)*(2*min.axis**2-2*min.axis*maj.axis+maj.axis**2)))
  x <- transform(x, a.hglass=(2*pi*(min.axis/2)**2))
  x <- transform(x, a.ellipse=(pi*min.axis*maj.axis))
  x <- transform(x, a.circle=(pi*(min.axis/2)**2))
  
  x <- transform(x, elongatedness=(maj.axis**2/min.axis))
  x <- transform(x, eccentricity=sqrt(1-((min.axis/2)/(maj.axis/2))**2))
  x <- transform(x, flattening=(maj.axis-min.axis)/(maj.axis+min.axis))
  x <- transform(x, perim.by.a=(perim/a.tot))
  x <- transform(x, perim.by.perimellipse=(perim/perim.ellipse))
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
  return(x)
}


#### get.fit.variables ####
## Extract columns for classification
get.fit.variables <- function(x, vars='m'){
  f.channels <- substr(unlist(regmatches(names(x), gregexpr(paste0("(f.tot.)+([a-z]{1}$)"), names(x)))), 7, 8)
  
  #id.vars <- c('pos', 'cellID', 'time', 'ucid', 't.frame')
  id.vars <- c('ucid','t.frame')
  
  #morph.vars <- c('xpos', 'ypos', 'fft.stat', 'perim', 'maj.axis', 'min.axis', 'rot.vol', 'con.vol', 'a.surf', 'sphere.vol')
  morph.vars <- c('fft.stat', 'perim', 'maj.axis', 'min.axis', 'a.surf', 'sphere.vol')
  
  #area.vars <- c('a.tot', 'a.tot.p1', 'a.tot.m1', 'a.tot.m2', 'a.tot.m3')
  area.vars <- c('a.tot')
  
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
  
  x <- make.fit.variables(x)
  
  new.vars <- names(x)
  new.vars <- new.vars[!(new.vars%in%default.vars)]
  
  all.vars <- names(x)
  
  ## Defines whether fluorescence ('f') and/or morphology ('m') variables should be extracted
  extr.vars <- c()
  if(sum(grep("m",vars))){
    morph.vars <- c(id.vars, morph.vars, area.vars, new.vars)
    extr.vars <- c(extr.vars,morph.vars)
  }
  if(sum(grep("f",vars))){
    fluor.vars <- c(f.vars, bg.vars, bg.f.vars)
    extr.vars <- c(extr.vars,fluor.vars)
  }
  #if(vars%in%c('fm', 'mf')){
  #  extr.vars <- c(id.vars, morph.vars, area.vars, f.vars, bg.vars, bg.f.vars, new.vars)
  #}
  #else if(vars=='m'){
  #  extr.vars <- c(id.vars, morph.vars, area.vars, new.vars)
  #}
  
  return(x[,names(x)%in%extr.vars])
}


#### k.means.classification ####
## Perform k-means classification to Cell-ID data. Appends a 'k' column indicating the assigned class, and a 'k.dist' column indicating the calculated Euclidean distance to the nearest centroid.
k.means.classification <- function(x, k=10, max.iterations=100, vars=NULL){
  ## Remove any pre-existing 'k' and 'k.dist' columns
  x$data <- x$data[,!(names(x$data)%in%c('k','k.dist'))]
  
  ## Get data to classify
  cellk <- subset(x$data, qc==TRUE)
  
  ## Create and get training variables
  cellk <- get.fit.variables(cellk, 'm')
  
  if(!is.null(vars)){
    vars <- c(c('ucid','t.frame','tag.type'),vars)
    cellk <- cellk[,names(cellk)%in%vars]
  }
  
  ## Remove rows with NAs
  n.cols <- ncol(cellk)
  na.idx <- unique(unlist(sapply(1:n.cols,function(x,i) which(is.na(x[,i])), x=cellk)))
  if(length(na.idx)) cellk <- cellk[-na.idx,]
  
  ## Extract ucids (id) and data (Xd)
  fit.vars <- names(cellk)[!(names(cellk)%in%c('ucid','t.frame','tag.type'))]
  id <- cellk[,c('ucid','t.frame')]
  Xd <- cellk[,fit.vars]
  
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
  
  while(i.diff>0 && i.count<max.iterations){
    ## Calculate distances of each row in Xd to each centroid in k.means
    dists <- apply(k.means,1,function(x,k) sqrt(rowSums(sweep(x, 2, k)**2)),x=Xd)
    
    ## Find nearest centroid for each row in Xd
    k.labs <- apply(dists,1,function(x) which(x==min(x))[1])
    i.diff <- sum(k.labs!=k.labs.old)
    i.count <- i.count + 1
    print(paste("Iteration:",i.count,"| Updates:",i.diff))
          
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