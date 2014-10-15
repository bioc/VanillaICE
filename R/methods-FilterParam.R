#' @param probability  minumum probability for the call
#' @param numberFeatures  minumum number of SNPs/nonpolymorphic features in a region
#' @param seqnames  the seqnames (character string or \code{Rle} to keep)
#' @param state  character:  the HMM states to keep
#' @param width the minimum widht of a region
#' @export
#' @rdname FilterParam-class
#' @seealso \code{\link{cnvFilter}} \code{\link{cnvSegs}} \code{\link{hmm2}}
#' @examples
#' ## To select CNV segments for which
#' ## - the CNV call has a 'posterior' probability of at least 0.95
#' ## - the number of features is at least 10
#' ## - the HMM states are 1 (homozygous deletion) or 2 (hemizygous deletion)
#' FilterParam(probability=0.95, numberFeatures=10, state=c("1", "2"))
#'
FilterParam <- function(probability=0.99, numberFeatures=10, seqnames=paste0("chr", c(1:22, "X", "Y")),
                        state=as.character(1:6), width=1L){
  new("FilterParam", probability=probability, numberFeatures=numberFeatures,
      seqnames=seqnames, state=state, width=1L)
}


#' @param object a \code{FilterParam} object
#' @aliases probability,FilterParam-method
#' @rdname FilterParam-class
setMethod("probability", "FilterParam", function(object) object@probability)

#' @aliases state,FilterParam-method
#' @rdname FilterParam-class
setMethod("state", "FilterParam", function(object) object@state)

setMethod("numberFeatures", "FilterParam", function(object) object@numberFeatures)

setMethod("seqnames", "FilterParam", function(x) x@seqnames)
setMethod("chromosome", "FilterParam", function(object) as.character(seqnames(object)))
setMethod("width", "FilterParam", function(x) x@width)

#' @aliases show,FilterParam-method
#' @rdname FilterParam-class
setMethod("show", "FilterParam", function(object){
  cat("An object of class 'FilterParam'\n")
  cat("   min. posterior probability of CNV call:", probability(object), "\n")
  cat("   min. no. of markers spanned by segment:", numberFeatures(object), "\n")
  cat("   min. width of segment                 :", width(object), "\n")
  states <- paste0(state(object), collapse=", ")
  cat("   selected HMM states                   :", states, "\n")
  chroms <- chromosome(object)
  if(length(chroms) > 6){
    chroms <- paste(paste(chroms[1:6], collapse=", "), "...")
  } else paste(chroms, collapse=", ")
  cat("   selected seqnames                     :", chroms, "\n")
})
