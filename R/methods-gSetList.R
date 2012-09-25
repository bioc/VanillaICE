
setMethod("[", signature(x="gSetList"),
	  function(x, i, j, ..., drop=FALSE){
		  ## using 'i' to subset markers does not really make
		  ## sense
		  ##
		  ## Use i to subset the list. example, x[1] is still a BeadStudioSetList, but is one chromosome
		  ##
		  if(!missing(i) & !missing(j)){
			  browser()
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  tmp <- lapply(tmp, function(x, j) {
					  x[, j, drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
			  x@phenoData <- phenoData(x)[j, ]
		  }
		  if(!missing(i) & missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
		  }
		  if(missing(i) & !missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- lapply(ad[[elt]], function(x, j) {
					  x[, j, drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@phenoData <- phenoData(x)[j, ]
		  }
		  return(x)
	  })



setMethod("show", signature(object="gSetList"),
	  function(object){
		  nm <- ls(assayData(object))[[1]]
		  lo <- length(assayData(object)[[nm]])
		  cat(class(object), " of length ", lo, "\n", sep="")
	  })

setMethod("assayDataList", signature(object="gSetList"),
	  function(object)  object@assayDataList)

setMethod("assayDataList", signature(object="oligoSetList"),
	  function(object)  object@assayDataList)

setMethod("featureDataList", signature(object="gSetList"),
	  function(object)  object@featureDataList)

setMethod("featureNames", signature(object="gSetList"),
	  function(object){
		  lapply(featureDataList(object), featureNames)
	  })

setMethod("position", signature(object="gSetList"),
	  function(object){
		  lapply(featureDataList(object), position)
	  })

setMethod("chromosome", signature(object="gSetList"),
	  function(object, as.list=FALSE, ...){
		  ##lapply(object, chromosome)
		  if(!as.list) object@chromosome else chromosomeList(object)
	  })

setMethod("varLabels", signature(object="gSetList"),
	  function(object) varLabels(phenoData(object)))

setMethod("pData", signature(object="gSetList"),
	  function(object) pData(phenoData(object)))

setMethod("$", signature(x="gSetList"),
	  function(x, name){
		  eval(substitute(phenoData(x)$NAME_ARG, list(NAME_ARG=name)))
	  })
setMethod("order", signature(...="gSetList"),
	  function(..., na.last=TRUE,decreasing=FALSE){
		  x <- list(...)[[1]]
		  for(i in seq_along(x)){
			  x[[i]] <- chromosomePositionOrder(x[[i]])
		  }
		  return(x)
	  })

##setMethod("genomeBuild", signature(object="gSetList"), function(object) object@genomeBuild)

setAs("gSetList", "list", function(from){
	to <- vector("list", length(from))
	for(i in seq_along(from)){
		to[[i]] <- from[[i]]
	}
	return(to)
})

setMethod("sapply", signature(X="gSetList"),
	  function(X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE){
		  listobject <- as(X, "list")
		  sapply(listobject, FUN, ..., simplify=simplify, USE.NAMES=USE.NAMES)
	  })

setMethod("dims", signature(object="gSetList"), function(object){
	nchr <- length(chromosome(object))
	nr <- sum(sapply(featureData(object), nrow))
	ns <- ncol(object[[1]])
	ds <- c(nr, nchr, ns)
	names(ds) <- c("features (total)", "list elements", "samples")
	return(ds)
})

setMethod("sampleNames", signature(object="gSetList"),
	  function(object) sampleNames(phenoData(object)))

setMethod("phenoData", signature(object="gSetList"),
	  function(object) object@phenoData)


setMethod("assayData", signature(object="gSetList"),
	  function(object) assayDataList(object))

setMethod("storageMode", "gSetList", function(object) storageMode(assayData(object)))

setMethod("length", signature(x="gSetList"), function(x) length(x@chromosome))

setReplaceMethod("$", "gSetList", function(x, name, value) {
	phenoData(x)[[name]] = value
	x
})

setReplaceMethod("phenoData",
                 signature=signature(
		 object="gSetList",
                   value="AnnotatedDataFrame"),
                 function(object, value) {
			 object@phenoData <- value
			 object
                 })

setReplaceMethod("pData",
                 signature=signature(
		 object="gSetList",
                   value="data.frame"),
                 function(object, value) {
			 pd <- phenoData(object)
			 pData(pd) <- value
			 phenoData(object) <- pd
			 object
                 })

setReplaceMethod("sampleNames", signature(object="gSetList", value="character"),
		 function(object,value){
			 pd <- phenoData(object)
			 sampleNames(pd) <- value
			 object@phenoData <- pd
			 return(object)
		 })


setReplaceMethod("featureData", signature(object="gSetList", value="list"),
		 function(object, value){
			 object@featureDataList <- value
			 object
		 })

setMethod("featureData", signature(object="gSetList"),
		 function(object){
			 object@featureDataList
		 })

setMethod("fvarLabels", signature(object="gSetList"),
	  function(object){
		  varLabels(featureData(object)[[1]])
	  })
