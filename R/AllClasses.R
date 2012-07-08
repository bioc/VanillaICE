setClass("HmmOptionList", contains="list")

setClass("Vit",
	 representation(emission="matrix",
			initialProb="numeric",
			transitionProb="matrix",
			arm="integer",
			numberStates="integer",
			numberFeatures="integer",
			viterbiStatePath="integer",
			normal2altered="numeric",
			altered2normal="numeric",
			altered2altered="numeric",
			normalIndex="integer",
			"VIRTUAL"))

setClass("Viterbi", contains="Vit",
	 representation(delta="matrix",
			pAA="numeric"))

setClass("Viterbi2", contains="Vit",
	 representation(forwardVariable="matrix",
			backwardVariable="matrix",
			scaleFactor="numeric"))



##	 representation(assayDataList="AssayData",
##			phenoData="AnnotatedDataFrame",
##			featureDataList="list",
##			chromosome="integer",
##			annotation="character",
##			genomeBuild="character"))
##setClass("BeadStudioSetList",
##	 representation(assayDataList="AssayData",
##			phenoData="AnnotatedDataFrame",
##			featureDataList="list",
##			chromosome="integer",
##			annotation="character",
##			genomeBuild="character"))
##

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




