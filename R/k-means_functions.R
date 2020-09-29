#' standardize_data
#'
#' Standardize data column-wise by calculating the respective Z-scores, i.e., subtracting the mean and dividing by the standard deviation.
#'
#' @param x a data.frame with numeric variables suitable for normalization
#'
#' @return a data.frame with normalized data
#'
standardize_data <- function(x){
  x <- t((t(x)-apply(x, 2, mean))/apply(x, 2, sd))
  rownames(x) <- NULL
  colnames(x) <- NULL
  cbind(rep(1,dim(x)[1]),x) # add intercept column
}

#' standardise_data
#'
#' Standardise data column-wise by calculating the respective Z-scores, i.e., subtracting the mean and dividing by the standard deviation.
#'
#' @param x a data.frame with numeric variables suitable for normalization
#'
#' @return a data.frame with normalized data
#'
standardise_data <- function(x){standardize_data(x)}


#' calculate_extra_morpho_vars
#'
#' Create custom fit variables
#'
#' @param x a cell.data data.frame
#'
#' @return a cell.data data.frame to which new variables have been added
#'
calculate_extra_morpho_vars <- function(x){
  x %>%
    mutate(perim.bud=(pi*maj.axis),
           perim.hglass=(pi*min.axis),
           perim.circle=(pi*(min.axis/2)**2)) %>% 
    mutate(a.bud=((pi/4)*(2*min.axis**2-2*min.axis*maj.axis+maj.axis**2)),
           a.hglass=(2*pi*(min.axis/2)**2),
           a.ellipse=(pi*min.axis*maj.axis),
           a.circle=(pi*(min.axis/2)**2)) %>% 
    mutate(elongatedness=(maj.axis**2/min.axis),
           eccentricity=sqrt(1-((min.axis/2)/(maj.axis/2))**2),
           flattening=(maj.axis-min.axis)/(maj.axis+min.axis),
           perim.by.a=(perim/a.tot),
           perim.by.ellipseperim=(perim/ellipse.perim),
           perim.by.perimbud=(perim/perim.bud),
           perim.by.spherevol=(perim/sphere.vol)) %>% 
    mutate(a.by.abud=(a.tot/a.bud),
           surf.by.vol=a.surf/sphere.vol,
           surf.by.perim=a.surf/perim,
           fit.var1=fft.stat/elongatedness,
           fit.var2=fft.stat/a.surf,
           fit.var3=fft.stat/a.tot,
           fit.var4=fft.stat/sphere.vol,
           fit.var5=a.surf/elongatedness,
           fit.var6=sphere.vol/elongatedness)
}


#' get_fit_vars
#'
#' Extract specific variables relevant for clustering
#'
#' @param x a cell.data data.frame
#' @param var_cats a character vector indicating whether morphology ('morpho') and/or fluorescence ('fluor') variables should be used for clustering.
#'
#' @return a reduced cell.data data.frame with selected clustering variables only
#'
get_fit_vars <- function(x, f.channels, var_cats=NULL, custom_vars=NULL){
  # Define ID variables; these are used to uniquely identify each row
  id.vars <- c('ucid','t.frame')
  
  ## Collect variables
  default.vars <- names(x)
  selected.vars <- c(id.vars,custom_vars)
  
  ## Conditionally include morphology variables
  if("morpho"%in%var_cats){
    # Set default morphology variables
    morpho.vars <- c('a.tot','fft.stat', 'perim', 'maj.axis', 'min.axis', 'a.surf', 'sphere.vol', 'ellipse.perim')
    
    # Calculate and append custom morphology variables
    x <- calculate_extra_morpho_vars(x)
    all.vars <- names(x)
    extra.morpho.vars <- all.vars[!(all.vars%in%default.vars)]
    
    selected.vars <- c(selected.vars,morpho.vars, extra.morpho.vars)
  }
  
  ## Conditionally include fluorescence variables
  if("fluor"%in%var_cats){
    # Select default fluorescence variables
    fluor.vars <- c("f.tot","f.bg","f.local.bg","f.nucl.tag1","f.nucl.tag3","f.nucl.tag5","f.vacuole","f.bg","f.tot.m1","f.tot.m3")
    # consider also adding "f"; this requires fixing Rcell2 to automatically calculate "f.t" and "cf.t" variables (if applicable).
    
    fluor.vars <- as.vector(sapply(fluor.vars,function(x,y) paste(x,y,sep="."),y=f.channels,simplify=TRUE))
    fluor.vars <- c(fluor.vars)
    selected.vars <- c(selected.vars,fluor.vars)
  }
  
  x %>% select(all_of(selected.vars))
  
  #return(x[,names(x) %in% selected.vars])
}

#' K-means clustering 
#'
#' Perform k-means clustering on a Cell-ID data.
#'
#' @param x cell.data object or a cell.data data.frame
#' @param k either a non-negative integer setting the number of clusters, or a data.frame with \code{ucid} and \code{t.frame} pairs giving the rows to be used as starting centroids.
#' @param max_iter non-negative integer. The maximum number of iterations allowed.
#' @param resume logical. If \code{TRUE} the algorithm picks up clustering from pre-assigned clusters found in a column with name \code{k} (default), or by the column name passed by the \code{tag_col} argument.
#' @param tag_col optional string for specifying the name of the column containing pre-defined clusters to be used when \code{resume} is set to \code{TRUE}. Defaults to \code{k}.
#' @param var_cats optional character vector specifying if morphological (\code{morpho}) and/or fluorescence (\code{fluor}) variables should be used for clustering.
#' @param custom_vars optional character vector specifying custom variables to be included for clustering.
#' 
#' @example 
#' 
#' @return a cell.data object with columns 'k' and 'k.dist' appended to the data.frame, indicating the assigned class and Euclidean distance to the assigned centroid, respectively.
#' @export
#'
kmeans_clustering <- function(x, k=10, max_iter=100, resume=FALSE, tag_col = 'k', var_cats=c("morpho"),custom_vars=NULL){
  ## Get filtered data
  if(is.cell.data(x)){
    xdata <- x$data[x$data$qc,]
    f.channels <- x$channels$posfix
    }
  else{
    xdata <- x
    f.channels <- substr(unlist(regmatches(names(xdata), gregexpr(paste0("(f.tot.)+([a-z]{1}$)"), names(xdata)))), 7, 8)
  }
  
  ## Extract columns for computing clustering
  cdata <- get_fit_vars(xdata,
                        f.channels,
                        var_cats,
                        custom_vars)
  
  ## Remove any rows with NAs
  na.idx <- unique(unlist(sapply(1:ncol(cdata),function(x,i) which(is.na(x[,i])), x=cdata)))
  if(length(na.idx)){
    cdata <- cdata[-na.idx,]
    xdata <- xdata[-na.idx,]
  }
  
  ## Extract id variables (id) and clustering data (cdata)
  id.vars <- c("ucid","t.frame")
  fit.vars <- names(cdata)[!(names(cdata)%in%id.vars)]
  id <- cdata %>% select(any_of(id.vars))
  cdata <- cdata %>% select(any_of(fit.vars))
  
  ## Standardize (normalize) clustering data
  cdata <- standardize_data(cdata)
  
  ### Run k-means classification ####
  if(resume){
    ### If clustering takes off from or resumes previous tagging
    ## Verify that tag_col exists
    stopifnot("missing k column" = tag_col%in%names(xdata))
    
    ## Extract tag_col
    tags <- xdata %>% pull(!!tag_col)
    unique_tags <- sort(unique(tags)[!is.na(unique(tags))])
    
    # Ensure numeric, consecutive cluster labels
    k.labs.current <- match(tags,unique_tags)
    k.labs.unique <- unique(k.labs.current[!is.na(k.labs.current)])
    
    ## Assign random cluster value to any missing values (NA)
    na.idx <- which(is.na(tags))
    k.labs.sampled <- sample(k.labs.unique,length(k.labs.current),replace=TRUE)
    k.labs.current[na.idx] <- k.labs.sampled[na.idx]
    
    ## Calculate current cluster centroids
    k.means <- t(sapply(1:length(unique_tags), function(i,x,n) colMeans(x[which(n==i),]), x=cdata, n=k.labs.current))
    
    ## Check for empty clusters, and randomly assign rows
    missing <- (1:k)[!(1:k%in%unique(k.labs.current))]
    k.sample <- sample(dim(cdata)[1],length(missing))
    k.means <- rbind(k.means,cdata[k.sample,])
  }
  else{
    ### Else, if tagging de novo
    ## Randomly pick k samples as initial cluster centroids
    k.sample <- sample(dim(cdata)[1],k)
    k.means <- cdata[k.sample,]
    
    ## Set random initial cluster assignments
    k.labs.current <- sample(1:k,dim(id)[1],replace=TRUE)
  }
  
  ## Loop over centroid calculation as long as labels continue getting updated
  ## Limit loop to max_iter iterations
  cat("Running k-means clustering...\n")
  for(i.count in 1:max_iter){
    ## Calculate distances of each row in cdata to each centroid in k.means
    dists <- apply(k.means,1,function(x,k) sqrt(rowSums(sweep(x, 2, k)**2)),x=cdata)
    
    ## Find nearest centroid for each row in cdata
    k.labs <- apply(dists,1,function(x) which(x==min(x))[1])
    i.diff <- sum(k.labs!=k.labs.current)
          
    ## Calculate new centroids
    k.means <- t(sapply(1:k, function(i,x,n) colMeans(x[which(n==i),]), x=cdata, n=k.labs))
    k.labs.current <- k.labs
    
    ## Print progress
    cat("\r  Iteration #",i.count," (",max_iter,")",". Updates: ",i.diff," ",sep="")
    flush.console() 
    
    if(i.diff==0){
      cat("\nk-means clustering converged successfully.\n")
      break
    }
    i.count <- i.count + 1
  }
  if(i.diff>0){
    cat("\nk-means clustering reached max number of iterations without converging.\nAppending clusters calculated on last iteration.\n")
  }
  
  ## Sum up data by id, cluster, and distance to cluster centroid
  k.iter <- seq_len(k)
  dists <- apply(k.means,1,function(x,k) sqrt(rowSums(sweep(x, 2, k)**2)),x=cdata)
  k.idx <- lapply(k.iter,function(i,l) which(l==i), l=k.labs)
  k.dists <- lapply(k.iter,function(d,k,i) d[k[[i]],i],d=dists,k=k.idx)
  k.ucids <- lapply(k.iter,function(u,k,i) u[k[[i]],], u=id, k=k.idx)
  k.data <- mapply(function(ucid,k,k.dist) cbind(ucid,k,k.dist),k.dist=k.dists,ucid=k.ucids,k=k.iter,SIMPLIFY=FALSE)
  k.data <- lapply(k.data,function(x) x[order(x[,3]),])
  k.data <- do.call(rbind.data.frame, k.data)
  
  ## Remove any pre-existing 'k' and 'k.dist' columns, and append new columns
  xdata <- xdata %>% mutate(k = NULL, k.dist = NULL) %>% left_join(k.data,by=c("t.frame","ucid"))
  
  ## Return original object with appended columns
  if(is.cell.data(x)) x$data <- xdata
  else x <- xdata
  return(x)
}