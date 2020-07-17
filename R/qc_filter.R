#public
#filters cells, modifies QC variable
#ToDo: warn for posible outdated variables
#' qc_filter
#'
#' @param X
#' @param filter
#' @param subset
#'
#' @return cell.data object
#' @export
#'
#' @examples
qc_filter <- function(X, filter, subset=NULL){

    filter = substitute(filter)
    subset = substitute(subset)

    #initializing QC if required
    if(is.null(X$data$QC)) X$data$QC = rep(TRUE, times = dim(X$data)[1])

    #saving the old filter for undo vector
    QC.last = X$data$QC

    attributes(QC.last) <- NULL
    QC.attr = attributes(X$data$QC)

    #updating the QC filter
    if(is.null(subset))
        X$data$QC = X$data$QC & eval(filter, X$data)
    else
        X$data$QC = X$data$QC & ( eval(filter,X$data) | !eval(subset, X$data) )

    #trating NAs as FALSE
    X$data$QC[is.na(X$data$QC)] <- FALSE

    #adding the information for undos as attributes of QC
    QC.attr.names = names(QC.attr)
    QC.history.names = names(X$QC.history)

    if(is.null(QC.attr.names))
        hNum = "QC0001"
    else
        hNum = paste("QC",
                     formatC(1 + as.numeric(substring(max(QC.history.names), 3, 6)),
                             width = 4,
                             flag = "0"),
                     sep="")

    attr(X$data$QC, hNum) <- QC.last

    for(i in QC.attr.names) {
        attr(X$data$QC,i) <- QC.attr[[i]]
    }

    if(length(QC.attr.names) >= 10) {
        attr(X$data$QC, min(QC.attr.names)) <- NULL
    }

    cer = sum(!X$data$QC) / length(X$data$QC)

    #adding call to QC.history
    tmp<-list()
    tmp[[hNum]] <- list(type = "filter",
                        filter = filter,
                        undo = NA,
                        cumulative.exclusion.ratio = cer)

    if(is.null(subset)) {
        tmp[[hNum]]$subset = NA
    } else {tmp[[hNum]]$subset = subset}

    X$QC.history <- c(X$QC.history,tmp)

    cat("cumulative row exclusion: ", round(100*cer,1), "%\n", sep="")

    return(X)
}

#*************************************************************************#
#public
#removes the last applied QC filter
#ToDo: allow multiple undos with second argument
#ToDo: warn for posible outdated variables
#' qc_undo
#'
#' @param X
#'
#' @return cell.data object
#' @export
#'
#' @examples
qc_undo <- function(X){
    #browser()
    if(is.null(X$data$QC) || length(X$QC.history) == 0) stop("No QC variable\n")
    if(is.null(attributes(X$data$QC))) stop("No more undos available\n")

    hNum = paste("QC", formatC(1 + as.numeric(substring(max(names(X$QC.history)),3,6)),width=4,flag="0"),sep="")
    QC.attr=attributes(X$data$QC)
    QC.restore=max(names(QC.attr))
    X$data$QC<-QC.attr[[QC.restore]]
    QC.attr[[QC.restore]]<-NULL
    for(i in names(QC.attr)) attr(X$data$QC,i)<-QC.attr[[i]]

    #adding call to QC.history
    tmp<-list()
    tmp[[hNum]]<-list(type="undo",undo=QC.restore)
    X$QC.history<-c(X$QC.history,tmp)
    X$QC.history[[QC.restore]]$undo<-TRUE

    QCr=X$QC.history[[QC.restore]]
    cat("undoing filter",deparse(QCr[["filter"]])
        ,ifelse(class(QCr[["subset"]])=="call",paste("( on",deparse(QCr[["subset"]]),")"),"")
        ,"\n")

    return(X)
}

#*************************************************************************#
#public
#resets the QC filter and undo history
#ToDo: allow subset?
#' qc_reset
#'
#' @param X
#'
#' @return cell.data object
#' @export
#'
#' @examples
qc_reset <- function(X){

    QC.reseted=c()

    if(is.null(X$data$QC) || length(X$QC.history)==0) hNum="QC0001"
    else {
        hNum=paste("QC",formatC(1+as.numeric(substring(max(names(X$QC.history)),3,6)),width=4,flag="0"),sep="")
        for(i in names(X$QC.history))
            if(is.na(X$QC.history[[i]]$undo)){
                X$QC.history[[i]]$undo=TRUE
                QC.reseted=c(QC.reseted,i)
            }
    }

    X$data$QC=rep(TRUE,times=dim(X$data)[1])

    #adding call to QC.history
    tmp<-list()
    tmp[[hNum]]<-list(type="reset",undo=QC.reseted)
    X$QC.history<-c(X$QC.history,tmp)
    cat("resetting all filters\n")

    return(X)
}
