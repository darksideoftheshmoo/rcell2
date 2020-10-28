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
#' Perform k-means clustering on Cell-ID data.
#' 
#' K-means clusters data by assigning each row to the nearest cluster based on its Euclidean distance to the center (centroid) of all clusters. After assigning all rows, centroid positions are updated by calculating the column means of all rows assigned to each cluster. Row assignment and centroid updates are performed iteratively until the algorithm converges, i.e., no rows are re-assigned after centroid positions have been updated.
#' 
#' The number of clusters is defined by the parameter \code{k}, and clustering can be either completely unsupervised (\code{k} is a number only setting the desired number of clusters), or semi-supervised (\code{k} is a data.frame of \code{ucid} and \code{t.frame} pairs defining which rows/cells to choose as starting centroids). If unsupervised, starting centroids are chosen randomly by sampling \code{k} rows from data. Semi-supervised clustering can also be achieved by indicating a column of pre-defined labels assigned to a subset of rows, which will then be used to calculate the positions of the starting centroids.
#' 
#' Note that this algorithm does not guarantee to find the optimum.
#'
#' @param x cell.data object or a cell.data data.frame
#' @param k either a non-negative integer setting the desired number of clusters, or a data.frame with \code{ucid} and \code{t.frame} pairs specifying the rows in \code{x} to be used as starting centroids.
#' @param max_iter The maximum number of iterations allowed.
#' @param resume logical. If \code{TRUE} the algorithm picks up clustering from pre-assigned clusters found in a column with name \code{k} (default), or by the column name passed by the \code{label_col} argument.
#' @param label_col optional string specifying the column containing pre-defined clusters used when \code{resume} is set to \code{TRUE}. This overrides the default column \code{k}.
#' @param var_cats optional character vector specifying whether pre-defined sets of morphological (\code{morpho}) and/or fluorescence (\code{fluor}) variables should be included for clustering. If no value is given and \code{custom_vars} is empty, this defaults to \code{morpho}.
#' @param custom_vars optional character vector specifying custom variables to be included for clustering. These are added to any variable sets specified by \code{var_cats}.
#' 
#' @example 
#' No example yet.
#' @return Depending on the data type provided by \code{x}, either a cell.data object or a cell.data data.frame with appended columns \code{k} and \code{k.dist}, indicating the assigned cluster and Euclidean distance to the cluster centroid, respectively.
#' @export
#'
kmeans_clustering <- function(x, k=10, max_iter=100, resume=FALSE, label_col = 'k', var_cats=NULL,custom_vars=NULL){
  max_iter <- as.integer(max_iter)
  stopifnot("invalid max_iter value" = is.integer(max_iter) & max_iter>0)
  
  if(is.null(custom_vars) & is.null(var_cats)) var_cats <- "morpho"
  x_has_qc <- "qc" %in% names(x)
  
  ## Get filtered data
  if(is.cell.data(x)){
    ## Do this if x is a cell.data object
    if(x_has_qc) xdata <- x$data[x$data[,"qc"],]
    f.channels <- x$channels$posfix
    }else{
    ## Do this if x is a data.frame
    if(x_has_qc) xdata <- x[x[,"qc"],]
    f.channels <- substr(unlist(regmatches(names(xdata), gregexpr(paste0("(f.tot.)+([a-z]{1}$)"), names(xdata)))), 7, 8)
    # f.channels <- sub(pattern = "f.tot.([a-z])", replacement = "\\1", 
    #                   x = grep(pattern = "f.tot.[a-z]",
    #                            x = names(cdata),
    #                            value = TRUE))
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
    ## Verify that label_col exists
    stopifnot("missing tag column" = label_col%in%names(xdata),
              "expected a numeric value for k for processing pre-clustered data" = !is.data.frame(k))
    
    ## Extract label_col
    tags <- xdata %>% pull(!!label_col)
    unique_tags <- sort(unique(tags)[!is.na(unique(tags))])
    
    # Ensure numeric, consecutive cluster labels
    k.labs.current <- match(tags,unique_tags)
    k.labs.unique <- unique(k.labs.current[!is.na(k.labs.current)])
    
    ## Assign random cluster value to any missing values (NA)
    na.idx <- is.na(tags)
    
    ## Calculate current cluster centroids
    k.means <- t(sapply(1:length(k.labs.unique), function(i,x,n) colMeans(x[which(n==i),]), x=cdata[!na.idx,], n=k.labs.current[!na.idx]))
    
    ## Randomly reassign rows if pre-assigned clusters > k
    if(k<length(k.labs.unique)){
      reassigned <- which(k.labs.current>k)
      k.labs.current[reassigned] <- sample(1:k,length(reassigned),replace=TRUE)
    }
    
    ## Randomly assign rows as starting centroids for empty clusters
    if(k>length(k.labs.unique)){
      n.missing <- length((1:k)[!(1:k%in%unique(k.labs.current))])
      k.sample <- sample(dim(cdata)[1],n.missing)
      k.means <- rbind(k.means,cdata[k.sample,])
    }
  }else{
    ### Else, if tagging de novo
    if(is.data.frame(k)){
      ## Extract rows chosen as initial cluster centroids
      k <- k %>% mutate(k_labs = 1)
      centroids.idx <- id %>% mutate(i_row = seq_along(t.frame)) %>%
        left_join(k,by=c("ucid","t.frame")) %>% 
        filter(!is.na(k_labs)) %>% select(ucid,t.frame,i_row)
      k.means <- cdata[centroids.idx$i_row,]
      k <- dim(k)[1]
    }else{
      ## Randomly pick k samples as initial cluster centroids
      k.sample <- sample(dim(cdata)[1],k)
      k.means <- cdata[k.sample,]
    }
  }
  ## Set random initial cluster assignments
  ## This does not affect assignments, only prevents later error
  k.labs.current <- sample(1:k,dim(id)[1],replace=TRUE)
  
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
    cat("\r  Iteration #",i.count," (",max_iter,")",". Re-assignments: ",i.diff,"           ",sep="")
    flush.console() 
    
    if(i.diff==0){
      cat("\nk-means clustering converged successfully.\n")
      break
    }
    i.count <- i.count + 1
  }
  if(i.diff>0){
    cat("\nk-means clustering reached max number of iterations without converging.\nAppending clusters from last iteration.\n")
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
  
  ## Remove any pre-existing 'k' and 'k.dist' columns from x, and append new columns
  if(is.cell.data(x)){
    x$data <- x$data %>% mutate(k = NULL, k.dist = NULL) %>% left_join(k.data,by=c("t.frame","ucid"))
  }else{
    x <- x %>% mutate(k = NULL, k.dist = NULL) %>% left_join(k.data,by=c("t.frame","ucid"))
  }
  return(x)
}
                   
