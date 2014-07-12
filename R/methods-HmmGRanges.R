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



setMethod("HmmGRanges", "missing", function(states, feature_starts, feature_chrom, loglik){
  new("HmmGRanges")
})

setMethod("HmmGRanges", "integer", function(states, feature_starts, feature_chrom, loglik){
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
  states <- sapply(split(states, segfactor), unique)
  hgr$state <- states
  metadata(hgr) <- list(loglik=loglik)
  names(hgr) <- NULL
  hgr
})

setMethod("HmmGRanges", "Rle", function(states, feature_starts, feature_chrom, loglik){
  HmmGRanges(as.integer(states), feature_starts, feature_chrom, loglik)
})

setMethod(state, "HmmGRanges", function(object) object$state)

setMethod(statei, "HmmGRanges", function(object) as.integer(object$state))

setMethod(statef, "HmmGRanges", function(object) factor(statei(object), levels=c(1,2,3,4,5,6)))

setMethod(loglik, "HmmGRanges", function(object) metadata(object)[["loglik"]])

setReplaceMethod("state", "HmmGRanges", function(object, value) {
  mcols(object)$state <- value
  object
})
