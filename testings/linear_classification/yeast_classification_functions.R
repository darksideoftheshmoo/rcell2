require(EBImage)

#### get.fit.variables ####
## Extract columns for classification
get.fit.variables <- function(x, vars){
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
  f.vars <- c(f.vars, 'f.total.c', 'f.total.r', 'f.total.y')
  
  #bg.vars <- c('a.local.bg', 'a.local', 'a.local2.bg', 'a.local2')
  bg.vars <- c()
  
  #bg.f.vars <- c('f.bg', 'f.local.bg', 'f.local2.bg')
  #bg.f.vars <- c('f.bg', 'f.local.bg')
  #bg.f.vars <- sapply(bg.f.vars, function(x, y) paste(x,y,sep='.'), y=f.channels)
  bg.f.vars <- c()
  
  #standard.vars <- c('pos', 't.frame', 'cellID', 'ucid', 'time', 'xpos', 'ypos', 'a.tot', 'fft.stat', 'perim', 'maj.axis', 'min.axis', 'rot.vol', 'con.vol', 'a.vacuole', 'a.tot.p1', 'a.tot.m1', 'a.tot.m2', 'a.tot.m3', 'a.local.bg', 'a.local', 'a.local2.bg', 'a.local2', 'a.surf', 'sphere.vol', 'QC', 'f.tot.c', 'f.nucl.c', 'a.nucl.c', 'f.vacuole.c', 'f.bg.c', 'f.tot.p1.c', 'f.tot.m1.c', 'f.tot.m2.c', 'f.tot.m3.c', 'xpos.nucl.c', 'ypos.nucl.c', 'f.nucl1.c', 'f.nucl.tag1.c', 'a.nucl1.c', 'f.nucl2.c', 'f.nucl.tag2.c', 'a.nucl2.c', 'f.nucl3.c', 'f.nucl.tag3.c', 'a.nucl3.c', 'f.nucl4.c', 'f.nucl.tag4.c', 'a.nucl4.c', 'f.nucl5.c', 'f.nucl.tag5.c', 'a.nucl5.c', 'f.nucl6.c', 'f.nucl.tag6.c', 'a.nucl6.c', 'f.local.bg.c', 'f.local2.bg.c', 'f.tot.r', 'f.nucl.r', 'a.nucl.r', 'f.vacuole.r', 'f.bg.r', 'f.tot.p1.r', 'f.tot.m1.r', 'f.tot.m2.r', 'f.tot.m3.r', 'xpos.nucl.r', 'ypos.nucl.r', 'f.nucl1.r', 'f.nucl.tag1.r', 'a.nucl1.r', 'f.nucl2.r', 'f.nucl.tag2.r', 'a.nucl2.r', 'f.nucl3.r', 'f.nucl.tag3.r', 'a.nucl3.r', 'f.nucl4.r', 'f.nucl.tag4.r', 'a.nucl4.r', 'f.nucl5.r', 'f.nucl.tag5.r', 'a.nucl5.r', 'f.nucl6.r', 'f.nucl.tag6.r', 'a.nucl6.r', 'f.local.bg.r', 'f.local2.bg.r', 'f.tot.t', 'f.nucl.t', 'a.nucl.t', 'f.vacuole.t', 'f.bg.t', 'f.tot.p1.t', 'f.tot.m1.t', 'f.tot.m2.t', 'f.tot.m3.t', 'xpos.nucl.t', 'ypos.nucl.t', 'f.nucl1.t', 'f.nucl.tag1.t', 'a.nucl1.t', 'f.nucl2.t', 'f.nucl.tag2.t', 'a.nucl2.t', 'f.nucl3.t', 'f.nucl.tag3.t', 'a.nucl3.t', 'f.nucl4.t', 'f.nucl.tag4.t', 'a.nucl4.t', 'f.nucl5.t', 'f.nucl.tag5.t', 'a.nucl5.t', 'f.nucl6.t', 'f.nucl.tag6.t', 'a.nucl6.t', 'f.local.bg.t', 'f.local2.bg.t', 'f.tot.y', 'f.nucl.y', 'a.nucl.y', 'f.vacuole.y', 'f.bg.y', 'f.tot.p1.y', 'f.tot.m1.y', 'f.tot.m2.y', 'f.tot.m3.y', 'xpos.nucl.y', 'ypos.nucl.y', 'f.nucl1.y', 'f.nucl.tag1.y', 'a.nucl1.y', 'f.nucl2.y', 'f.nucl.tag2.y', 'a.nucl2.y', 'f.nucl3.y', 'f.nucl.tag3.y', 'a.nucl3.y', 'f.nucl4.y', 'f.nucl.tag4.y', 'a.nucl4.y', 'f.nucl5.y', 'f.nucl.tag5.y', 'a.nucl5.y', 'f.nucl6.y', 'f.nucl.tag6.y', 'a.nucl6.y', 'f.local.bg.y', 'f.local2.bg.y')
  
  default.vars <- names(x)
  
  x <- make.fit.variables(x)
  
  new.vars <- names(x)
  new.vars <- new.vars[!(new.vars%in%default.vars)]
  
  custom.vars <- c('tag.type','h1','h2','h3','h4','h5','h6','h7','h8')
  
  all.vars <- names(x)
  
  if(vars%in%c('fm', 'mf')){
    extr.vars <- c(id.vars, morph.vars, area.vars, f.vars, bg.vars, bg.f.vars, custom.vars, new.vars)
  }
  else if(vars=='m'){
    extr.vars <- c(id.vars, morph.vars, area.vars, custom.vars, new.vars)
  }
  
  return(x[,names(x)%in%extr.vars])
  
  #return(list(data=x, vars=extr.vars))
}


#### sigmoid ####
## Define sigmoid function
# Multiplication by .9999 to avoid ones
sigmoid <- function(x) (1/(1+exp(-x)))*.99999999

#### cost.function ####
## Define cost function
# With regularization
cost.function <- function(X, y, th, lambda){
  m <- dim(X)[1]
  z <- t(th%*%t(X))
  h <- sigmoid(z)
  return((1/m)*sum(y*log(h)-(1-y)*log(1-h)) + (lambda/(2*m))*sum(th**2))
}

#### gradients ####
## Calculate gradient (partial derivatives)
# With regularization
gradients <- function(X, y, th, lambda){
  m <- dim(X)[1]
  z <- t(th%*%t(X))
  h <- sigmoid(z)
  grad <- (1/m)*t(h-y)%*%X
  reg.grad <- (1/m)*t(h-y)%*%X + (lambda/m)*th
  grad[,-1] <- reg.grad[,-1]
  return(grad)
}

#### make.onehot ####
## Convert labels to onehot vectors
make.onehot <- function(y){
  n <- max(y)
  y_onehot <- c()
  for(i in y){
    onehot <- rep(0,n)
    onehot[i] <- 1
    y_onehot <- rbind(y_onehot, onehot)
  }
  rownames(y_onehot) <- NULL
  return(y_onehot)
}

#### X.normalize ####
## Normalize X
X.normalize <- function(X){
  X <- t((t(X)-apply(X, 2, mean))/apply(X, 2, sd))
  rownames(X) <- NULL
  colnames(X) <- NULL
  return(X)
}

#### gradient.descent ####
## Gradient descent
# Returns list with one element corresponding to a unique alpha/lambda combination.
gradient.descent <- function(X, y, alpha, lambda, iter){
  grad.results <- list()
  i_count <- 1
  for(a in alpha){
    for(l in lambda){
      i <- max(dim(y)[2],1)
      j <- dim(X)[2]
      th <- matrix(0, i, j)
      J_iter <- c()
      for(i_iter in 1:iter){
        J <- cost.function(X, y, th, l)
        grad <- gradients(X, y, th, l)
        th <- th - a*grad
        J_iter <- c(J_iter, J)
      }
      grad.results[[i_count]] <- list(J = J, theta = th, alpha = a, lambda = l, J_iter = J_iter)
      i_count <- i_count+1
    }
  }
  return(grad.results)
}

#### make.fit.variables ####
#Create custom fit variables
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

#### hu.moments ####
# Returns log Hu moments for a provided object.
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
  # The last Hu moment is mentioned in Wikipedia
  I8 <- eta(xy,1,1) * ((eta(xy,3,0) + eta(xy,1,2))**2 - (eta(xy,0,3) + eta(xy,2,1))**2) -
    (eta(xy,2,0) - eta(xy,0,2)) * (eta(xy,3,0) + eta(xy,1,2)) * (eta(xy,0,3) + eta(xy,2,1))
  
  #return(c(I1, I2, I3, I4, I5, I6, I7, I8))
  #hu <- c(I1, I2, I3, I4, I5, I6, I7, I8)
  hu <- abs(c(I1, I2, I3, I4, I5, I6, I7, I8))
  return(-sign(hu)*log(hu))
}

#### get.mask.data ####
# Analyzes objects in provided BF.out image, returns a data frame with variables
# calculated from masks, bound to their corresponding ucid.
# Takes a Cell-ID data frame 'x' and path to a BF.out image 'img.name'.
get.mask.data <- function(x, img.name){
  # Get image title for printing process
  img.title <- unlist(regmatches(img.name,gregexpr("([A-Z]{2,3})+(_Position)[[:digit:]]+(_{0,1}[a-z0-9]{0,15}.tif)", img.name)))
  print(paste0(img.title,'...'))
  
  min.obj.size <- 20
  bf.id <- get.pos.time(img.name)
  
  img <- imageData(create.mask(img.name))
  img.lab <- img
  img.lab[img!=.4] <- 0
  img.lab[img==.4] <- 1
  img.lab <- bwlabel(img.lab)
  
  kern = makeBrush(3, shape='box')
  obj.count <- max(img.lab)
  #obj.count <- 10 #Uncomment to limit number of objects analyzed
  
  #mask.data <- list()
  mask.data <- vector("list", obj.count)
  obj.rm <- c()
  for(i in 1:obj.count){
    #print(paste0(i,'/',obj.count))
    if(sum(img.lab==i)>=min.obj.size){
      img.single <- img
      img.dil <- img.lab==i
      img.dil <- dilate(img.dil, kern)
      img.single[img.dil==0] <- 0
      
      #contour.pixels <- which(img.single==1, arr.ind = TRUE)
      #colnames(contour.pixels) <- c('x','y')
      #fill.pixels <- which(img.single==.4, arr.ind = TRUE)
      #colnames(fill.pixels) <- c('x','y')
      cell.pixels <- which(img.single>0, arr.ind = TRUE)
      colnames(cell.pixels) <- c('x','y')
      xlim <- range(cell.pixels[,1])
      ylim <- range(cell.pixels[,2])
      xypos <- round(c(median(xlim),median(ylim)))
      hu <- hu.moments(cell.pixels)
      
      ##Pair with correct ucid
      pos.col <- match('pos',names(x))
      tf.col <- match('t.frame',names(x))
      sub.col <- match(c('ucid','xpos','ypos'),names(x))
      pos.idx <- which(x[,pos.col]==bf.id[1])
      tf.idx <- which(x[,tf.col]==bf.id[2])
      subset.idx <- intersect(pos.idx, tf.idx)
      xy <- x[subset.idx,sub.col]
      
      dists <- sqrt(rowSums(sweep(xy[,2:3], 2, xypos)**2))
      min.dist <- min(dists)
      ucid <- xy[which(dists==min.dist),1][1]
      
      mask.data[[i]] <- list(img = img.name,
                             pos = bf.id[1],
                             t.frame = bf.id[2],
                             ucid = ucid,
                             #tucid = as.numeric(paste0(ucid,bf.id[2])),
                             #cont.px = contour.pixels,
                             #fill.px = fill.pixels,
                             cell.px = cell.pixels,
                             #xmin = xlim[1],
                             #xmax = xlim[2],
                             #ymin = ylim[1],
                             #ymax = ylim[2],
                             #area = nrow(cell.pixels),
                             #xy.dist = min.dist,
                             h1 = hu[1],
                             h2 = hu[2],
                             h3 = hu[3],
                             h4 = hu[4],
                             h5 = hu[5],
                             h6 = hu[6],
                             h7 = abs(hu[7]),
                             h8 = abs(hu[8])
                             #xypos = xypos)
      )
    }
    else{
      #If 
      obj.rm <- c(obj.rm, i)
      mask.data[[i]] <- 0
    }
  }
  mask.data[obj.rm] <- NULL
  
  ##Find and remove excessively small objects
  #obj.rm <- c()
  #obj.size <- unlist(lapply(mask.data, function(x) nrow(x$fill.px)))
  #obj.rm <- which(obj.size<=5)
  #mask.data[obj.rm] <- NULL
  
  ##Find duplicated ucids
  obj.ucid <- unlist(lapply(mask.data, function(x) x[['ucid']]))
  obj.dist <- unlist(lapply(mask.data, function(x) x[['xy.dist']]))
  dup.ucid <- unique(obj.ucid[duplicated(obj.ucid)])
  obj.rm <- c()
  #For every duplicated ucid, keep only object closest to Cell-ID xypos
  for(i in dup.ucid){
    ucid.idx <- which(obj.ucid==i)
    dist.idx <- obj.dist[ucid.idx]
    obj.rm <- c(obj.rm,ucid.idx[which(dist.idx!=min(dist.idx))])
  }
  obj.rm <- unique(obj.rm)
  mask.data[obj.rm] <- NULL
  
  names(mask.data) <- unlist(lapply(mask.data, function(x) x$tucid))
  
  return(mask.data)
}

#### make.mask.df ####
## Condense mask.data list into data frame, discarding single pixel data
make.mask.df <- function(xlist){
  excl.vars <- c('img','cont.px','fill.px','cell.px')
  xlist.vars <- names(xlist[[1]])
  xlist.vars <- xlist.vars[-which(xlist.vars%in%excl.vars)]
  df <- as.data.frame(matrix(unlist(lapply(xlist, function(x) x[xlist.vars])), ncol=length(xlist.vars), byrow = TRUE))
  names(df) <- xlist.vars
  return(df)
}

#### get.img.names ####
## Loads image filenames, order by channel
get.img.names<-function(path){
  filenames <- as.data.frame(dir(path))
  img.type <- c('bf','tfp','yfp','rfp','cfp')
  img.type <- toupper(img.type)
  img.names <- c()
  for(i in img.type){
    #Extract filenames matching XX_PositionYYY_timeZZZ.tif
    img <- unlist(regmatches(filenames[,1], gregexpr(paste0("(",i,"_Position)[[:digit:]]+(_time)[[:digit:]]+(\\.tif$)"), filenames[,1])))
    img <- paste0(path,sort(img))
    
    #Extract filenames matching XX_PositionYYY_timeZZZ.tif.out.tif
    img.out <- unlist(regmatches(filenames[,1], gregexpr(paste0("(",i,"_Position)[[:digit:]]+(_time)[[:digit:]]+(\\.tif\\.out.tif$)"), filenames[,1])))
    img.out <- paste0(path,sort(img.out))
    
    if(length(img)!=length(img.out)){
      print(paste0(i,' and ',i,'.out file counts not matching.'))
    }
    if(length(img)){
      img.comb <- cbind(img, img.out)
      colnames(img.comb) <- c(i, paste0(i,'.out'))
      img.names <- cbind(img.names, img.comb)
    }
  }
  
  return(img.names)
}

#### get.pos.time ####
## Extracts position and timepoint for all BF.out filenames provided
get.pos.time <- function(x){
  #Extract positions from RegEx pattern
  ps <- as.integer(sub("_Position","",unlist(regmatches(x, gregexpr("(_Position)[[:digit:]]+", x)))))
  #Extract timepoints from RegEx pattern
  tm <- as.integer(sub("_time","",unlist(regmatches(x, gregexpr("(_time)[[:digit:]]+", x)))))
  
  x <- cbind(ps, tm-1)
  colnames(x) <- c('pos', 't.frame')
  
  return(x)
}

#### create.mask ####
## Fill all holes in a BF.out image.
#Fill color is different (0.4) from contour color (1.0)
create.mask <- function(img.name){
  img <- EBImage::readImage(img.name)
  img[img==max(img)] <- 1
  img[img!=1] <- 0
  img.contour <- img
  img <- fillHull(img)
  img[img==1] <- .4
  img[img.contour==1] <- 1
  return(img)
}

#### assign.tags.positions ####
# This function assigns tags to positions at a specific t.frame. The tag table
# should correspond to an image stack of one or more positions at the same t.frame.

# Takes a Cell-ID data frame 'cdata', a vector 'tfile' with tagging filepath(s),
# a vector 'p' with position(s) to which tags should be assigned, and optionally a
# value 'tm' indicating the t.frame used for assigning tags. If 'tm' is not provided,
# the default value t.frame==0 is assumed.
assign.tags.positions <- function(cdata, tfile, p, tm){
  if(missing(tm)) tm <- 0
  
  # Import tag table
  tdata<-read.table(tfile,sep=",",skip=1)[,1:4]
  
  # Convert tag table to matrix
  tdata <- as.matrix(tdata)
  
  # Get relevant Cell-ID variable indices
  cvars <- c('pos','t.frame','ucid','xpos','ypos')
  cvar.idx <- sapply(cvars,function(x,y) which(y==x), colnames(cdata))
  
  # Extract Cell-ID position data at t.frame==0, then separate into list of positions
  cdata <- cdata[intersect(which(cdata[,cvar.idx[1]]%in%p),
                           which(cdata[,cvar.idx[2]]==tm)),cvar.idx]
  
  #tf <- sort(unique(cdata[,2]))
  cdata <- lapply(p, function(x,y) y[which(y[,1]==x),], y=cdata)
  
  # Calculate distances between each tag and cell coordinates
  dists <- apply(tdata, 1, function(x,y){
    sqrt(rowSums(sweep(y[[x[2]]][,4:5], 2, x[3:4])**2))
  },
  y=cdata)
  
  # Extract minimum distance for each tag
  min.dists <- lapply(dists, min)
  
  # Get index for minimum distances
  min.dists.idx <- mapply(function(x,y) which(y==x), x=dists, y=min.dists)
  
  # Check for redundancies in min.dists.idx, keep only first
  min.dists.idx <- unlist(lapply(min.dists.idx, function(x){
    if(length(x)>1){
      print(x)
      x <- x[1]}
    x
  }
  ))
  
  # Get ucids
  closest.ucid <- sapply(seq_len(dim(tdata)[1]), function(x,y,z,d){
    y[[z[x]]][d[x],3]
  },
  y=cdata, z=tdata[,2], d=min.dists.idx)
  ###closest.ucid <- sapply(min.dists.idx, function(x,y) y[x,3], y=cdata)
  ###y[y[,2]==(x[2]-1),4:5]
  
  # Bind ucids to tag table and reorder
  tdata <- as.data.frame(cbind(tdata[,1:2],closest.ucid))
  tdata <- transform(tdata,t.frame=tm)[,c(4,3,1)]
  
  ###names(tdata) <- c('tag.type','ucid')
  names(tdata) <- c('t.frame','ucid','tag.type')
  
  return(tdata)
}
