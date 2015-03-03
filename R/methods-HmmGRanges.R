setValidity("HmmGRanges", function(object){
  nms <- colnames(mcols(object))
  if(!"numberFeatures" %in% nms) {
    return("'numberFeatures' not a variable in mcols()")
  }
  if(!"state" %in% nms) {
    return("'state' not a variable in mcols()")
  }
  TRUE
})


# @aliases HmmGRanges,missing-method
# @rdname HmmGRanges-class
setMethod("HmmGRanges", "missing", function(states, feature_starts,
                                            feature_chrom,
                                            loglik,
                                            emission_param=EmissionParam()){
  new("HmmGRanges")
})

setMethod("emissionParam", "HmmGRanges", function(object) object@emission_param)

setReplaceMethod("emissionParam", c("HmmGRanges", "EmissionParam"), function(object, value){
  object@emission_param <- value
  object
})

# @aliases HmmGRanges,integer-method
# @rdname HmmGRanges-class
setMethod("HmmGRanges", "integer", function(states, feature_starts, feature_chrom, loglik, emission_param=EmissionParam()){
  ## we can't split by states because the same state could span multiple chromosomes
  segfactor <- paste(feature_chrom, states, sep=";")
  segfactor <- Rle(segfactor)
  segment_index <- Rle(seq_along(runLength(segfactor)), runLength(segfactor))
  segfactor <- factor(as.integer(segment_index))
  feature_start_list <- split(feature_starts, segfactor)
  starts <- sapply(feature_start_list, min)
  ends <- sapply(feature_start_list, max)
  chr <- sapply(split(feature_chrom, segfactor), unique)
  if(!is.character(chr)) stop("one segment has multiple chromosomes, resulting in a list of chromosomes instead of a character vector.")
  gr <- GRanges(chr, IRanges(starts, ends),
                numberFeatures=elementLengths(feature_start_list))
  hgr <- as(gr, "HmmGRanges")
  emissionParam(hgr) <- emission_param
  states <- sapply(split(states, segfactor), unique)
  hgr$state <- states
  metadata(hgr) <- list(loglik=loglik)
  names(hgr) <- NULL
  hgr
})

# @aliases HmmGRanges,Rle-method
# @rdname HmmGRanges-class
setMethod("HmmGRanges", "Rle", function(states, feature_starts, feature_chrom, loglik){
  HmmGRanges(as.integer(states), feature_starts, feature_chrom, loglik)
})

#' Accessor for copy number state
#'
#' Extract the copy number state for each genomic interval.
#' @param object a \code{HmmGRanges} object
#' @aliases state,HmmGRanges-method
#' @rdname HmmGRanges-methods
setMethod(state, "HmmGRanges", function(object) object$state)

setMethod(statei, "HmmGRanges", function(object) as.integer(object$state))

setMethod(statef, "HmmGRanges", function(object) factor(statei(object), levels=c(1,2,3,4,5,6)))

setMethod(loglik, "HmmGRanges", function(object) metadata(object)[["loglik"]])

setReplaceMethod("state", "HmmGRanges", function(object, value) {
  mcols(object)$state <- value
  object
})

#' @aliases  cnvSegs,HmmGRanges-method
#' @rdname cnvFilter
setMethod("cnvSegs", "HmmGRanges", function(object, filters=FilterParam(state=as.character(c(1,2,5,6)))){
  cnv <- cnvFilter(object, filters)
  cnv
})

setMethod("numberFeatures", "HmmGRanges", function(object) object$numberFeatures)



.apply_cnv_filters <- function(g, filters){
  if(length(g)==0) return(g)
  keep <- width(g) > width(filters)
  keep <- keep & state(g) %in% state(filters)
  keep <- keep & numberFeatures(g) >= numberFeatures(filters)
  keep <- keep & chromosome(g) %in% seqnames(filters)
  keep <- keep & g$prCall >= probability(filters)
  g[keep]
}


setMethod("cnvFilter", "GRanges", function(object, filters=FilterParam()){
  .apply_cnv_filters(object, filters)
})

grangesData <- function(granges, se, expandFUN, ylim){
  size <- expandFUN(granges)
  id <- unique(granges$id)
  if(length(id) > 1) stop("only one sample can be passed to grangesData")
  chr <- as.character(seqnames(granges))
  starts <- as.integer(pmax(start(granges)-size, 1))
  ends <- as.integer(pmin(end(granges)+size, seqlengths(granges)[chr]))
  gexpanded <- GRanges(seqnames(granges), IRanges(starts, ends))
  hits <- findOverlaps(gexpanded, rowRanges(se))
  region <- queryHits(hits)
  i <- subjectHits(hits)
  j <- match(granges$id[1], colnames(se))
  ##df2 <- data.frame(y=c(lrr(se)[i, 1], baf(se)[i, 1]),
  df2 <- data.frame(y=c(lrr(se)[i, j], baf(se)[i, j]),
                    x=rep(start(se)[i], 2),
                    type=rep(c("LRR","BAF"), each=length(i)),
                    region=rep(region, 2),
                    state=rep(state(granges)[region], 2))
  df2$y <- threshold(df2$y, ylim)
  df2$starts <- start(granges)[region]
  df2$ends <- end(granges)[region]
  dflist <- split(df2, df2$region)
  dflist
}

.xyplotlist <- function(granges, se, param){
  ylimits <- yLimits(param)
  expandfun <- expandFun(param)
  dflist <- grangesData(granges, se, expandfun, ylim=ylimits[[2]])
  figlist <- lapply(dflist, function(dat){
    at <- pretty(range(dat$x), n=7)
    labels <- prettyNum(at/1000, big.mark=",")
    xyplot(y~x | type, dat, pch=20, col="gray40", cex=0.5,
           layout=c(1,2),
           scales=list(y=list(relation="free",rot=0), x=list(at=at, labels=labels, cex=0.7, axs="i")),
           panel=function(x, y, starts, ends, type, state, ..., subscripts){
             start <- starts[subscripts]
             end <- ends[subscripts]
             type <- (type[subscripts])[1]
             ##if(type=="LRR")  y <- 2*2^y
             panel.xyplot(x,y,...)
             panel.abline(v=c(start, end), lty=2)
             panel.abline(h=0)
           }, start=dat$starts, end=dat$ends,
           ylim=ylimits,
           type=dat$type,
           state=dat$state,
           xlab="Position (kb)", ylab="",
           strip=FALSE,
           strip.left=strip.custom(style=1, horizontal=FALSE))
  })
  figlist
}

#' @aliases xyplotList,HmmGRanges,SnpArrayExperiment-method
#' @rdname plotting
setMethod("xyplotList", c("HmmGRanges", "SnpArrayExperiment"), function(granges, se, param=HmmTrellisParam()){
  .xyplotlist(granges, se, param)
})

#' @aliases xyplotList,GRangesList,SnpArrayExperiment-method
#' @rdname plotting
setMethod("xyplotList", c("GRangesList", "SnpArrayExperiment"), function(granges, se, param=HmmTrellisParam()){
  se <- se[, names(granges)]
  figList <- lapply(granges,
                    function(cnv, se, param){
                      .xyplotlist(cnv, se, param)
                    }, se=se,
                    param=param)
  names(figList) <- colnames(se)
  figList
})

setMethod("reduce", "HmmGRanges", function(x, ...){
  g <- as(x, "GRanges")
  gr <- reduce(g, ...)
  hits <- findOverlaps(gr, g)
  nf <- split(numberFeatures(x), queryHits(hits))
  nf <- sapply(nf, sum)
  st <- split(state(x), queryHits(hits))
  st <- sapply(st, function(x) paste(unique(x), collapse=","))
  prCall <- sapply(split(x$prCall, queryHits(hits)), min)
  id <- sapply(split(x$id, queryHits(hits)), unique)
  gr$numberFeatures <- nf
  gr$state <- st
  gr$prCall <- prCall
  gr$id <- id
  as(gr, "HmmGRanges")
})
