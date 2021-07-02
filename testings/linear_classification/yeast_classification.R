#### SET PATHS #############################################

## Load functions
source("/Users/andreas/Library/Mobile Documents/com~apple~CloudDocs/Analys mikroskopi/functions_rcell_analysis.R")
source("/Users/andreas/Library/Mobile Documents/com~apple~CloudDocs/Analys mikroskopi/yeast_classification_functions.R")

## Set experiment date and title
set.paths(
  expfolder="190528_memoria_wt_preinduction_wo_inhibitor"
  #,inroot="/Volumes/HFS+/Andreas" #Custom image file location
  #,inroot="/Users/andreas/Mikroskopi" #Custom image file location
  #,outroot="/Users/andreas/Moln/pCloud/Jobb/Buenos Aires/Experiment pCloud" #Custom output location
)->expdata

D<-load.filtered(expdata)


## Erase previous tagging and filtering ####

D$data <- transform(D$data,filter=QC)
D$data$QC <- T
tag.col <- names(D$data)%in%'tag.type'
D$data <- D$data[,!tag.col]


#### GET MASK DATA #########################################

## Get BF.out filenames and select image to mask
img.names <- get.img.names(expdata$inpath)

## Select images to analyze
#pos: 1:4, 11:20, 29:32
img.sel <- c(1:4,11:20,29:32)*7-6
#img.sel <- c(1,2)

## Get and merge mask variables
mask.df <- lapply(img.sel,function(i,x,nm){
  bf <- nm[i,'BF.out']
  m.list <- get.mask.data(x, bf)
  m.df <- make.mask.df(m.list)
  return(m.df)
  #x <- merge(x, m.df, all.x=TRUE)
},x=D$data,nm=img.names)
mask.df <- do.call(rbind.data.frame, mask.df)
D$data <- merge(D$data, mask.df, all.x=TRUE)


#### TAGGING ###############################################

#1 Shmoo, 2 Broad projection, 3 Arrested, 4 G1, 5 S, 6 G2/M,
#7 Cell junk, 8 Background

tag.path <- paste0(expdata$inpath,'ML Stacks/')
tag.table <- 'CellCounter_BF_class_stack_results.csv'
tag.file <- paste0(tag.path,tag.table,sep='')
tag.positions <- c(11:20,29:30)
ucid.tags <- assign.tags.positions(D$data, tag.file, tag.positions)
D$data <- merge(D$data, ucid.tags, all.x=TRUE)


#### EXTRACT DATA FOR TRAINING ##########################################

# Get data to classify
untagged.idx <- which(D$data[,'t.frame']==0 & D$data[,'pos']%in%c(1:4,11:20,29:32))
tagged.idx <- which(!is.na(D$data[,'tag.type']))
subset.idx <- unique(c(untagged.idx,tagged.idx))
celld <- D$data[subset.idx,]

# Create and get training variables
celld <- get.fit.variables(celld, 'm')

# Give NAs in tag.type class value 0
celld[is.na(celld[,'tag.type']),'tag.type'] <- 0

# Remove any rows with NAs
n.cols <- ncol(celld)
na.idx <- unique(unlist(sapply(1:n.cols,function(x,i) which(is.na(x[,i])), x=celld)))
celld <- celld[-na.idx,]

## Remove duplicated rows
celld <- celld[!duplicated(celld[,'ucid']),]

# Adjust classes for training, save into new variable 'tag', and remove
celld <- transform(celld, tag=0)
celld[celld$tag.type==1,'tag'] <- 1 #shmoo
celld[celld$tag.type==3,'tag'] <- 2 #Arrested
celld[celld$tag.type==4,'tag'] <- 3 #G1
celld[celld$tag.type==5,'tag'] <- 4 #S
celld[celld$tag.type==6,'tag'] <- 5 #G2/M
celld[celld$tag.type==8,'tag'] <- 6 #background

celld <- celld[,!(names(celld)=='tag.type')]
celld <- transform(celld, tag.type=tag)
celld <- celld[,!(names(celld)=='tag')]

# Show cells from each class
i.tag = 6
sum(celld$tag.type==i.tag)
show.ucid <- subset(celld, tag.type==i.tag)$ucid
#cimage(D,~cell,channel="BF.out",N=49,subset=ucid%in%show.ucid & t.frame==0,box.size=40)


#### SETUP TRAINING ####################################

## Extract ucids, data (X), and labels (y)
fit.vars <- names(celld)[!(names(celld)%in%c('ucid','tag.type'))]
id <- celld[,'ucid']
X <- celld[,fit.vars]
y <- celld[,'tag.type']

## Normalize X
X <- X.normalize(X)

do.pca <- 0
if(do.pca){
  # Principle Component Analysis
  pca <- svd(var(X))
  U <- pca$u
  Dval <- pca$d
  V <- pca$v
  
  # Project data, and show contribution of each component to total variance
  Z <- X%*%U
  var.tot <- sum(diag(var(X%*%U)))
  var.cont <- Dval/var.tot
  round(cumsum(var.cont),3)
  qplot(x=1:length(var.cont),y=cumsum(var.cont),ylim=0:1)
  n.dim <- which(cumsum(var.cont)>=.999)[1]
  
  # Project data
  U.red <- U[,1:n.dim]
  Z <- X%*%U.red
  
  # Recover dimensionality reduced data
  #X <- Z%*%t(U.red)
  
  X <- Z
}

## Add intercept term
X <- cbind(rep(1, nrow(X)), X)

## Make onehot vectors. Unclassified cells (tag.type==0) get all-zero elements
y <- make.onehot(y)

## Make balanced training and validation sets
# Class sizes are equalized to that of the smallest class
tags <- seq_len(dim(y)[2])
tag.idx <- lapply(tags,function(i,y) which(y[,i]==1), y=y)
tag.size <- sapply(tag.idx,length)
tag.sampling <- sapply(tag.idx,function(x,n) sample(x,n), n=min(tag.size))
tag.train <- sample(as.vector(tag.sampling[1:round(min(tag.size)*.8),]))
tag.valid <- sample(as.vector(tag.sampling[-(1:round(min(tag.size)*.8)),]))
  
X.train <- X[tag.train,]
X.valid <- X[tag.valid,]
X.test <- X[-as.vector(tag.sampling),]
y.train <- y[tag.train,]
y.valid <- y[tag.valid,]
y.test <- y[-as.vector(tag.sampling),]
id.train <- id[tag.train]
id.valid <- id[tag.valid]
id.test <- id[-as.vector(tag.sampling)]


#### RUN TRAINING ####################################################

## Run gradient descent
alpha <- 1/(2**(0:3))
alpha <- .1
lambda <- 1/(2**(0:3))
lambda <- 1/(2**4)
iter <- 1000

grad.results <- gradient.descent(X.train, y.train, alpha, lambda, iter)

## Extract results from gradient descent
grad.results.cost <- unlist(lapply(grad.results, function(x) x$J))
cost.results <- lapply(grad.results, function(x) x$J_iter)
cost.iter <- t(ldply(cost.results, .fun=function(x, y) rbind(y, x), y=c()))
cost.iter <- melt(cost.iter)
names(cost.iter) <- c('iteration', 'group', 'cost')
qplot(x=iteration, y=cost, color=factor(group), data=cost.iter)

# Predict on validation set
th <- grad.results[[1]]$th
pred.valid <- sigmoid(X.valid%*%t(th))
pred.valid <- pred.valid/rowSums(pred.valid) # Convert to probabilities
rownames(pred.valid) <- id.valid
pred.valid.amax <- apply(pred.valid, 1, function(x) which(x==max(x)))

rownames(y.valid) <- id.valid
y.valid.amax <- apply(y.valid, 1, function(x) which(x==max(x)))

# Show accuracy
sum(pred.valid.amax==y.valid.amax)/dim(X.valid)[1]

# Show prediction vs ground truth by ucid
pred.valid.summary <- cbind(pred.valid.amax,y.valid.amax)

# Show incorrect predictions
i.tag <- 2
incorrect.pred <- id.valid[t(diff(t(pred.valid.summary)))!=0]
predicted.class <- id.valid[pred.valid.amax==i.tag]
show.ucid <- incorrect.pred[incorrect.pred%in%predicted.class]
#cimage(D,~cell,channel="BF.out",N=49,subset=ucid%in%show.ucid & t.frame==0,box.size=40)

# Show prediction scores for incorrectly predicted cells
pred.scores <- cbind(round(pred.valid,3), pred.valid.amax, y.valid.amax)
colnames(pred.scores) <- c(1:dim(y.valid)[2], 'Pred','GT')
pred.scores[id.valid%in%show.ucid,]

# Show correct predictions
i.tag <- 3
correct.pred <- id.valid[t(diff(t(pred.valid.summary)))==0]
predicted.class <- id.valid[pred.valid.amax==i.tag]
show.ucid <- correct.pred[correct.pred%in%predicted.class]
#cimage(D,~cell,channel="BF.out",N=49,subset=ucid%in%show.ucid & t.frame==0,box.size=40)

# Show prediction scores for correctly predicted cells
pred.scores <- cbind(round(pred.valid,3), pred.valid.amax, y.valid.amax)
colnames(pred.scores) <- c(1:dim(y.valid)[2], 'Pred','GT')
pred.scores[id.valid%in%show.ucid,]

# Plot weights
th.named <- th[,-1]
fit.vars <- names(celld)
colnames(th.named) <- fit.vars[!(fit.vars%in%c('ucid','tag.type'))]
rownames(th.named) <- 1:dim(th.named)[1]
th.named <- melt(th.named)
names(th.named) <- c('group', 'variable', 'value')

f1<-function(x) return(0)
qplot(x=variable, y=value, color=factor(group), data=th.named, ylim=c(-3,3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_function(fun=f1, color="red", linetype='dashed', size=.3)

# Predict on test set (untagged)
pred.test <- sigmoid(X.test%*%t(th))
pred.test <- pred.test/rowSums(pred.test) # Convert to probabilities
rownames(pred.test) <- id.test
pred.test.amax <- apply(pred.test, 1, function(x) which(x==max(x)))

# Show predictions
i.tag <- 6
show.ucid <- id.test[pred.test.amax==i.tag]
cimage(D,~cell,channel="BF.out",N=25,subset=ucid%in%show.ucid & t.frame==0,box.size=25)

## Train on complete (balanced) tagged set
grad.results <- gradient.descent(rbind(X.train,X.valid), rbind(y.train,y.valid), alpha, lambda, iter)
th <- grad.results[[1]]$th

# Predict on test set (untagged)
pred.test <- sigmoid(X.test%*%t(th))
pred.test <- pred.test/rowSums(pred.test) # Convert to probabilities
rownames(pred.test) <- id.test
pred.test.amax <- apply(pred.test, 1, function(x) which(x==max(x)))

# Show predictions
i.tag <- 1
show.ucid <- id.test[pred.test.amax==i.tag]
cimage(D,~cell,channel="BF.out",N=25,subset=ucid%in%show.ucid & t.frame==0,box.size=25)


#### ASSIGN PREDICTIONS #################################################

# Make data frame with predicted tags
pred.ucid.tags <- cbind(id.test, pred.test.amax)
rownames(pred.ucid.tags) <- NULL
colnames(pred.ucid.tags) <- c('ucid','tag')
pred.ucid.tags <- as.data.frame(pred.ucid.tags)

# Make data frame with GT tags
gt.ucid.tags <- celld[celld[,'tag.type']>0,names(celld)%in%c('ucid','tag.type')]
colnames(gt.ucid.tags) <- c('ucid','tag')

# Erase overlapping ucids from predicted tags
pred.ucid.tags <- pred.ucid.tags[!pred.ucid.tags[,1]%in%gt.ucid.tags[,1],]

# Combine GT and predicted tags into single data frame
ucid.tags <- rbind(pred.ucid.tags,gt.ucid.tags)

# Merge to new variable 'tag' in cell data
celld <- merge(celld, ucid.tags, all.x = TRUE)

# Set any untagged cells to 'tag==0'
celld[celld[,'tag']==0,'tag'] <- 0


#### K-MEANS ###########################################################

## Get data to classify
subset.idx <- which(D$data[,'t.frame']==0)
cellk <- D$data[subset.idx,!names(D$data)=='tag.type']
#cellk <- D$data[subset.idx,]
#cellk[is.na(cellk[,'tag.type']),'tag.type'] <- 0

## Create and get training variables
cellk <- get.fit.variables(cellk, 'm')

## Remove any rows with NAs
n.cols <- ncol(cellk)
na.idx <- unique(unlist(sapply(1:n.cols,function(x,i) which(is.na(x[,i])), x=cellk)))
if(length(na.idx)) cellk <- cellk[-na.idx,]

## Extract ucids (id) and data (X)
fit.vars <- names(cellk)[!(names(cellk)%in%c('ucid','tag.type'))]

#tags <- cellk[,'tag.type']
id <- cellk[,'ucid']
X <- cellk[,fit.vars]


#### Run k-means ###################################################
k <- 10

## Normalize X
X <- X.normalize(X)

## Randomly pick k samples as initial centroids
k.sample <- sample(dim(X)[1],k)
k.means <- X[k.sample,]

## (Optional) Calculate initial centroids from 'tag.type' means
k.tags <- sort(unique(tags))[-1]
k <- length(k.tags)
k.means <- t(sapply(k.tags,function(i,x,tags) colMeans(x[tags==i,]), x=X, tags=tags))

## Set random initial centroid assignments
k.labs.old <- sample(1:k,length(id),replace=TRUE)

lab.diff <- length(k.labs.old)
while(lab.diff>0){
  ## Calculate minimum distances of each row in X to k means
  dists <- apply(k.means,1,function(x,k) sqrt(rowSums(sweep(x, 2, k)**2)),x=X)
  
  ## Find closest centroid for each row in X
  k.labs <- apply(dists,1,function(x) which(x==min(x))[1])
  lab.diff <- sum(k.labs!=k.labs.old)
  print(lab.diff)
  
  ## Calculate new centroids
  k.means <- t(sapply(1:k, function(i,x,n) colMeans(x[which(n==i),]), x=X, n=k.labs))
  k.labs.old <- k.labs
}

## Show class sizes
sapply(sort(unique(k.labs)), function(k,l) sum(l==k), l=k.labs)

## Sum up data by ucid, centroid, and distance to centroid
k.iter <- seq_len(k)
dists <- apply(k.means,1,function(x,k) sqrt(rowSums(sweep(x, 2, k)**2)),x=X)
k.idx <- lapply(k.iter,function(i,l) which(l==i), l=k.labs)
k.dists <- lapply(k.iter,function(d,k,i) d[k[[i]],i],d=dists,k=k.idx)
k.ucids <- lapply(k.iter,function(u,k,i) u[k[[i]]], u=id, k=k.idx)
k.data <- mapply(function(ucid,k,dist) cbind(ucid,k,dist),dist=k.dists,ucid=k.ucids,k=k.iter)
k.data <- lapply(k.data,function(x) x[order(x[,3]),])
k.data <- do.call(rbind.data.frame, k.data)

## Plot distribution of distances to centroid for each class
qplot(x=dist,data=k.data,facets=~k,binwidth=.1,fill=factor(k))

## Show cells from specific class
i.k <- 8
qplot(x=dist,data=subset(k.data,k==i.k),facets=~k,binwidth=.1,fill=factor(k),xlim=c(0,20))
show.ucid <- k.data[k.data$k==i.k,'ucid']
cimage(D,~cell,channel="BF.out",N=49,subset=ucid%in%show.ucid & t.frame==0,box.size=25)

## Show cells within class with smallest and largest distances to centroid
k.subset <- subset(k.data,k==i.k)
show.ucid <- k.subset[1:12,1]
cimage(D,~cell,channel="BF.out",N=49,subset=ucid%in%show.ucid & t.frame==0,box.size=25)
show.ucid <- k.subset[order(-k.subset[,3])[1:12],1]
cimage(D,~cell,channel="BF.out",N=49,subset=ucid%in%show.ucid & t.frame==0,box.size=25)

## Collect data to be erased
k.erase <- subset(k.data,ucid%in%show.ucid)
cellk <- subset(cellk,!(ucid%in%k.erase$ucid))

## Select class for re-analysis
k.select <- subset(k.data,k==7)
k.idx <- which(cellk$ucid%in%k.select$ucid)
id <- cellk[k.idx,'ucid']
X <- cellk[k.idx,fit.vars]


#### PCA #########################################

## Run PCA on each centroid class and project down onto two dimensions
# The centroid in each class will be at (0,0), and the objects assigned
# to that centroid scattered in 2D along their two principal components
X.red <- list()
id.red <- list()
k.iter <- seq_len(k)
for(i.k in k.iter){
  ## Get data belonging to a certain centroid
  k.select <- subset(k.data,k==i.k)
  k.idx <- which(cellk$ucid%in%k.select$ucid)
  id <- cellk[k.idx,'ucid']
  X <- cellk[k.idx,fit.vars]
  
  X <- X.normalize(X)
  
  # Principle Component Analysis
  pca <- svd(var(X))
  U <- pca$u
  Dval <- pca$d
  V <- pca$v
  
  # Project data, and show contribution of each component to total variance
  Z <- X%*%U
  var.tot <- sum(diag(var(X%*%U)))
  var.cont <- Dval/var.tot
  round(cumsum(var.cont),3)
  qplot(x=1:length(var.cont),y=cumsum(var.cont),ylim=0:1)
  n.dim <- which(cumsum(var.cont)>=.999)[1]
  
  # Project data
  #U.red <- U[,1:n.dim]
  U.red <- U[,1:2]
  Z <- X%*%U.red
  
  # Recover dimensionality reduced data
  #X <- Z%*%t(U.red)
  
  X.red[[i.k]] <- as.data.frame(cbind(id,Z))
}

X.red <- mapply(function(x,i) transform(x,k=i), x=X.red, i=k.iter, SIMPLIFY = FALSE)
X.red <- do.call(rbind.data.frame, X.red)
names(X.red) <- c('ucid','x','y','k')

qplot(y=y,x=x,data=X.red,facets=~k,color=factor(k),size=I(.8))

## Extract and display some cells of interest
show.ucid <- subset(X.red, k==4 & y < -10)$ucid
cimage(D,~cell,channel="BF.out",N=49,subset=ucid%in%show.ucid & t.frame==0,box.size=25)