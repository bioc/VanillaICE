# Constructor for HMM class
#
# The contructor \code{HMM} creates and object of class \code{HMM}.


#' @param granges a \code{GRanges} object
#' @param param  a \code{HmmParam} object
#' @param posterior matrix of posterior probabilities
#' @param filters  an object of class \code{FilterParam}
#' @export
#' @rdname HMM
HMM <- function(granges=GRanges(), param=HmmParam(), posterior=matrix(), filters=FilterParam()){
  new("HMM", granges=granges, param=param, posterior=posterior, filters=filters)
}

setMethod("granges", "HMM", function(x, use.mcols=FALSE, ...){
  x@granges
})

setMethod("filters", "HMM", function(object) object@filters)


setMethod("posterior", "HMM", function(object) object@posterior)

setMethod("getHmmParams", "HMM", function(object) object@param)
setMethod("emissionParam", "HMM", function(object) emissionParam(object@param))

setReplaceMethod("granges", "HMM", function(x, value){
  x@granges <- value
  x
})

setReplaceMethod("seqinfo", "HMM",
                 function(x, value){
                   seqinfo(granges(x)) <- value
                   x
                 })



setReplaceMethod("emission", "HMM", function(object, value){
  hmm_param <- object@param
  emission(hmm_param) <- value
  object@param <- hmm_param
  object
})

#' @aliases cnvSegs,HMM-method
#' @rdname cnvFilter
setMethod("cnvSegs", "HMM", function(object, filters=FilterParam(state=as.character(c(1,2,5,6)))){
  x <- segs(object)
  cnv <- cnvFilter(x, filters)
  cnv
})

setMethod("cnvFilter", "HMM", function(object, filters=FilterParam()){
  granges <- cnvSegs(object)
  .apply_cnv_filters(granges, filters)
})

#' @aliases state,HMM-method
#' @rdname HMM
setMethod("state", "HMM", function(object) segs(object)$state)



setMethod("numberFeatures", "HMM", function(object) segs(object)$numberFeatures)

setMethod("duplication", "HMM", function(object, filters=FilterParam(state=c("5","6"))){
  cnvFilter(object, filters)
})

setMethod("hemizygous", "HMM", function(object, filters=FilterParam(state="2")){
  cnvSegs(object, filters)
})

setMethod("homozygous", "HMM", function(object, filters=FilterParam(state="1")){
  cnvSegs(object, filters)
})

setMethod("deletion", "HMM", function(object, filters=FilterParam(state=c("1","2"))){
  cnvSegs(object, filters)
})

#' @param object a \code{HMM} object
#' @aliases show,HMM-method
#' @rdname HMM
#' @export
setMethod("show", "HMM", function(object){
  cat("Object of class 'HMM'\n")
  cat("  no. segments:", length(granges(object)), "\n")
  dups <- duplication(object)
  hemi <- hemizygous(object)
  homo <- homozygous(object)
  cat("  no. duplications        :", length(dups), "\n")
  cat("  no. hemizygous deletions:", length(hemi), "\n")
  cat("  no. homozygous deletions:", length(homo), "\n")
  cat("See posterior(), segs(), cnvSegs(), getHmmParams(), \n")
})

setMethod("segs", "HMM", function(object) object@granges)
