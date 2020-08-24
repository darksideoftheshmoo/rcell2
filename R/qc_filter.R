#public

#' qc_filter
#'
#' Filters cells that DO NOT MATCH the condition in \code{filter}. i.e. if you want
#' to REMOVE cells with a.tot < 100, call \code{qc_filter(X, a.tot > 100)}.
#' In addition, subsets selects a subset of cells to which apply the filter.
#' i.e. to apply the same area filter but only to cells with f.tot.y values under 50,
#' call \code{qc_filter(X, a.tot > 100, subset = f.tot.y < 50)}
#'
#' It does not delete the filtered rows, it only modifies the qc variable. To effectively
#' remove them, use \link{qc_execute}.
#'
#' @param X cell.data object
#' @param filter conditions to filter
#' @param subset subset of cells to select for filtering
#'
#' @return cell.data object
#' @export
#'
#' @examples
qc_filter <- function(X, filter, subset=NULL){

    filter = substitute(filter)
    subset = substitute(subset)

    #initializing qc if required
    if(is.null(X$data$qc)) X$data$qc = rep(TRUE, times = dim(X$data)[1])

    #saving the old filter for undo vector
    qc.last = X$data$qc

    attributes(qc.last) <- NULL
    qc.attr = attributes(X$data$qc)

    #updating the qc filter
    if(is.null(subset))
        X$data$qc = X$data$qc & eval(filter, X$data)
    else
        X$data$qc = X$data$qc & ( eval(filter,X$data) | !eval(subset, X$data) )

    #trating NAs as FALSE
    X$data$qc[is.na(X$data$qc)] <- FALSE

    #adding the information for undos as attributes of qc
    qc.attr.names = names(qc.attr)
    qc.history.names = names(X$qc.history)

    if(is.null(qc.attr.names))
        hNum = "qc0001"
    else
        hNum = paste("qc",
                     formatC(1 + as.numeric(substring(max(qc.history.names), 3, 6)),
                             width = 4,
                             flag = "0"),
                     sep="")

    attr(X$data$qc, hNum) <- qc.last

    for(i in qc.attr.names) {
        attr(X$data$qc,i) <- qc.attr[[i]]
    }

    if(length(qc.attr.names) >= 10) {
        attr(X$data$qc, min(qc.attr.names)) <- NULL
    }

    cer = sum(!X$data$qc) / length(X$data$qc)

    #adding call to qc.history
    tmp<-list()
    tmp[[hNum]] <- list(type = "filter",
                        filter = filter,
                        undo = NA,
                        cumulative.exclusion.ratio = cer)

    if(is.null(subset)) {
        tmp[[hNum]]$subset = NA
    } else {tmp[[hNum]]$subset = subset}

    X$qc.history <- c(X$qc.history,tmp)

    cat("cumulative row exclusion: ", round(100*cer,1), "%\n", sep="")

    return(X)
}

#*************************************************************************#
#public

#' qc_undo
#'
#' Undoes the last qc_filter applied to the cell.data object.
#'
#' @param X cell.data
#'
#' @return cell.data
#' @export
#'
#' @examples
qc_undo <- function(X){
    #browser()
    if(is.null(X$data$qc) || length(X$qc.history) == 0) stop("No qc variable\n")
    if(is.null(attributes(X$data$qc))) stop("No more undos available\n")

    hNum = paste("qc", formatC(1 + as.numeric(substring(max(names(X$qc.history)),3,6)),width=4,flag="0"),sep="")
    qc.attr=attributes(X$data$qc)
    qc.restore=max(names(qc.attr))
    X$data$qc<-qc.attr[[qc.restore]]
    qc.attr[[qc.restore]]<-NULL
    for(i in names(qc.attr)) attr(X$data$qc,i)<-qc.attr[[i]]

    #adding call to qc.history
    tmp<-list()
    tmp[[hNum]]<-list(type="undo",undo=qc.restore)
    X$qc.history<-c(X$qc.history,tmp)
    X$qc.history[[qc.restore]]$undo<-TRUE

    qcr=X$qc.history[[qc.restore]]
    cat("undoing filter",deparse(qcr[["filter"]])
        ,ifelse(class(qcr[["subset"]])=="call",paste("( on",deparse(qcr[["subset"]]),")"),"")
        ,"\n")

    return(X)
}

#*************************************************************************#
#public
#resets the qc filter and undo history
#ToDo: allow subset?
#' qc_reset
#'
#' resets the \code{qc} variable to \code{TRUE}, removing all applied \code{qc_filters}.
#'
#' @param X cell.data object
#'
#' @return cell.data object
#' @export
#'
#' @examples
qc_reset <- function(X){

    qc.reseted=c()

    if(is.null(X$data$qc) || length(X$qc.history)==0) hNum="qc0001"
    else {
        hNum=paste("qc",formatC(1+as.numeric(substring(max(names(X$qc.history)),3,6)),width=4,flag="0"),sep="")
        for(i in names(X$qc.history))
            if(is.na(X$qc.history[[i]]$undo)){
                X$qc.history[[i]]$undo=TRUE
                qc.reseted=c(qc.reseted,i)
            }
    }

    X$data$qc=rep(TRUE,times=dim(X$data)[1])

    #adding call to qc.history
    tmp<-list()
    tmp[[hNum]]<-list(type="reset",undo=qc.reseted)
    X$qc.history<-c(X$qc.history,tmp)
    cat("resetting all filters\n")

    return(X)
}


#*************************************************************************#

#ToDo: modify undo value of prevoius qc.history elements
#ToDo: use subset to code this function
if(getRversion() >= "2.15.1") utils::globalVariables(c("qc"))
#' qc_execute
#'
#' removes cells with \code{qc == False}.
#'
#' @param X cell.data
#'
#' @return cell.data
#' @export
#'
#' @examples
qc_execute <- function(X){
    qc.attr=attributes(X$data$qc)
    qc.attr.names=names(qc.attr)
    qc.history.names=names(X$qc.history)
    if(is.null(qc.attr.names))
        hNum="qc0001"
    else
        hNum=paste("qc",formatC(1+as.numeric(substring(max(qc.history.names),3,6)),width=4,flag="0"),sep="")

    #calculating the cummulative row exclusion before deleting the registers
    cer=sum(!X$data$qc)/length(X$data$qc)

    cat("Eliminating ",format(round(100*cer,1),digits=3,nsmall=1),"% of the dataset registers\n",sep="")

    X$data<-subset(X$data,qc)
    X$data$qc=rep(TRUE,times=dim(X$data)[1])

    tmp<-list()
    tmp[[hNum]]<-list(type="execute",filter=NA,undo=FALSE,cumulative.exclusion.ratio=cer,subset=NA)

    #setting previous filters as definitive
    X$qc.history<-
        lapply(X$qc.history,FUN=function(l){
            if(is.na(l$undo)) l$undo<-FALSE
            return(l)
        })

    X$qc.history<-c(X$qc.history,tmp)

    return(X)
}
