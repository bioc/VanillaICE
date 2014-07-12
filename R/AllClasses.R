setClass("HmmOptionList", contains="list")


setClass("Viterbi",
         representation(state="integer",
                        loglik="numeric",
                        forward_backward="matrix"))

setValidity("BeadStudioSetList", function(object){
	nms <- ls(assayData(object))
	if(!all(c("baf", "lrr") %in% nms)){
		msg <- "baf and lrr are required elements of the assayData"
		return(msg)
	}
	if(length(object) > 0){
		msg <- validAssayDataDims(assayData(object))
		if(!all(msg == TRUE)) return(msg)
		elt <- (ls(assayDataList(object)))[[1]]
		b <- assayDataList(object)[[elt]]
		if(length(chromosome(object)) != length(b)){
			return("chromosome slot must be the same length as the length of the list for each assayData element")
		}
	}
	if(!identical(sampleNames(object), sampleNames(phenoData(object)))){
		stop("sampleNames of 'BeadStudioSetList' object must be the same as the sampleNames of the phenoData")
	}
	if(length(featureDataList(object)) != length(chromosome(object))){
		return("each chromosome should have an element in the featureDataList")
	}
	if(length(featureDataList(object)) > 0){
		featureDataClasses <- sapply(featureDataList(object), class)
		if(!unique(featureDataClasses) == "GenomeAnnotatedDataFrame"){
			return("featureDataList must be comprised of GenomeAnnotatedDataFrame(s)")
		}
	}
})

setClass("EmissionParam",
         representation=representation(
           cn_means="numeric",
           cn_sds="numeric",
           baf_means="numeric",
           baf_sds="numeric",
           initial="numeric",
           EMupdates="integer",
           CN_range="numeric",
           proportionOutlier="numeric",
           temper="numeric"))

setClass("TransitionParam",
         representation(
           taup="numeric",
           taumax="numeric"))

## should change the name to ViterbiLogLik
setClass("LogLik", representation(loglik="numeric", tolerance="numeric"))

setClass("HmmParam",
         representation(
           emission="matrix",
           emission_param="EmissionParam",
           transition="numeric",
           chromosome="character",
           loglik="LogLik",
           viterbi="Viterbi",
           verbose="logical"))


#' @export
setClass("SnpDataFrame", contains="DataFrame")

#' @export
setClass("SnpGRanges", contains="GRanges",
         representation(elementMetadata="SnpDataFrame"))


#' @rdname SnpArrayExperiment-class
#' @export
setClass("SnpArrayExperiment", contains = "SummarizedExperiment",
         representation(rowData="SnpGRanges"))

#' @aliases HmmGRanges-class
#' @rdname HmmGRanges-class
#' @export
setClass("HmmGRanges", contains="GRanges")

##setClass("SummarizedHMM",
##         representation(loglik="numeric",
##                        rowData="HmmGRanges"))
##

#' @export
setClass("ArrayViews",
         representation(colData="DataFrame",
                        rowData="GRanges",
                        index="integer",
                        filePaths="character",
                        scale="numeric",
                        datadir="character"))

#' @export
setClass("GStudioViews", contains="ArrayViews")

#' @export
setClass("GStudioScanParams", representation(readfun="function",
                                             index_genome="integer",
                                             cnvar="character",
                                             bafvar="character",
                                             gtvar="character",
                                             scale="numeric",
                                             select="integer"))
