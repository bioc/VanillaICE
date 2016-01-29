setValidity("HMMList", function(object){
  msg <- TRUE
  if(length(object) > 0){
    cls <- sapply(object, class)
    all_elements_HMM <- all(cls=="HMM")
    if(!all_elements_HMM) msg <- "All elements must be an HMM object"
  }
  msg
})

#' Constructor for \code{HMMList} class
#'
#' The constructor function for the \code{HMMList} class. The
#' constructor is useful for representing a list of \code{HMM}
#' objects.
#'
#' @param object a list. Each element of the list is in instance of
#' the \code{HMM} class.
#' @seealso \code{\linkS4class{HMMList}} \code{\linkS4class{HMM}} \code{\link{hmm2}}
#' @export
#' @rdname HMMList
HMMList <- function(object){
  if(missing(object))   return(new("HMMList") )
  if(is(object, "list")){
    return(new("HMMList", object))
  }
}

#' @param object a \code{HMMList} object
#' @export
#' @rdname HMMList-class
setMethod("show", "HMMList", function(object){
  cat("An object of class 'HMMList'\n")
  cat("  Length: ", length(object), "\n")
  cat("  Use [[j]] to retrieve results for the jth sample. \n")
  cat("  Use unlist to coerce to a 'GRanges' object. \n")
})

#' @param x a \code{HMMList} object
#' @param recursive logical;  currently ignored
#' @param use.names logical;  currently ignored
#' @rdname HMMList-class
setMethod("unlist", "HMMList", function(x, recursive=TRUE, use.names=TRUE){
  grl <- lapply(x, granges)
  id <- names(x)
  if(!is.null(id)){
    nm <- NULL
    grl <- foreach(g=grl, nm=id) %do% {
      g$id <- nm
      g
    }
  }
  g <- unlist(GRangesList(grl))
  g
})

#' @aliases segs,HMMList-method
#' @rdname cnvFilter
setMethod("segs", "HMMList", function(object) unlist(object))

#' @aliases hemizygous,HMMList-method
#' @rdname cnvFilter
setMethod("hemizygous", "HMMList", function(object) {
  x <- GRangesList(lapply(object, hemizygous))
  id <- rep(names(x), elementNROWS(x))
  x <- unlist(x)
  x$id <- id
  x
})

#' @aliases homozygous,HMMList-method
#' @rdname cnvFilter
setMethod("homozygous", "HMMList", function(object) {
  x <- GRangesList(lapply(object, homozygous))
  id <- rep(names(x), elementNROWS(x))
  x <- unlist(x)
  x$id <- id
  x
})

#' @aliases duplication,HMMList-method
#' @rdname cnvFilter
setMethod("duplication", "HMMList", function(object) {
  x <- GRangesList(lapply(object, duplication))
  id <- rep(names(x), elementNROWS(x))
  x <- unlist(x)
  x$id <- id
  x
})

#' @aliases cnvSegs,HMMList-method
#' @rdname cnvFilter
setMethod("cnvSegs", "HMMList", function(object, filters=FilterParam(state=as.character(c(1,2,5,6)))){
  x <- unlist(object)
  cnvFilter(x, filters)
})

#' @aliases cnvFilter,HMMList-method
#' @rdname cnvFilter
setMethod("cnvFilter", "HMMList", function(object, filters=FilterParam()){
  x <- unlist(object)
  cnvFilter(x, filters)
})

computeHeterozygosity <- function(hmmlist, snpexp, filter.param=FilterParam()){
  if(ncol(snpexp) != length(hmmlist)) stop("Number of elements in the hmmlist must be the same as the number of columns in snpexp ")
  j <- NULL
  results <- foreach(x = hmmlist, j=seq_along(hmmlist)) %do% {
    g <- cnvFilter(x, filter.param)
    se2 <- snpexp[, j]
    if(length(g)==0 ){
      g$percent.het <- numeric()
      granges(x) <- g
      return(x)
    }
    g$percent.het <- rep(NA, length(g))
    hits <- findOverlaps(g, se2)
    hitlist <- split(hits, queryHits(hits))
    for(i in seq_along(hitlist)){
      h <- hitlist[[i]]
      k <- subjectHits(h)
      m <- unique(queryHits(h))
      g$percent.het[m] <- round(mean(isHeterozygous(se2[k, ], cutoff=c(0.4, 0.6))), 2)
    }
    granges(x) <- g
    return(x)
  }
  HMMList(results)
}

flagHemDelWithHets <- function(segs, thr=0.01){
  flag <- rep(0L, length(segs))
  is_hemizygous <- state(segs) %in% "2"
  is_het <- segs$percent.het >= thr
  flag[is_hemizygous & is_het] <- 1L
  flag
}

flagHomDelWithHets <- function(segs, thr=0.2){
  flag <- rep(0L, length(segs))
  is_hom <- state(segs) %in% "1"
  is_het <- segs$percent.het >= thr
  flag[is_hom & is_het] <- 1L
  flag
}
