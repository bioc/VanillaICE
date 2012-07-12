setMethod("hmm", signature(object="oligoSnpSet"),
	  function(object, ...){
		  if("baf" %in% assayDataElementNames(object)){
			  return(hmmBeadStudioSet(object, ...))
		  } else{
			  return(hmmOligoSnpSet(object, ...))
		  }
	  })

setMethod("hmm", signature(object="BeadStudioSet"),
	  function(object, ...){
		  hmmBeadStudioSet(object, ...)
	  })

setMethod("hmm", signature(object="SnpSet2"),
	  function(object, ...){
		  hmmSnpSet(object, ...)
	  })

setMethod("hmm", signature(object="BeadStudioSetList"),
	  function(object, ...){
		  hmmBeadStudioSetList(object, ...)
	  })

setMethod("hmm", signature(object="oligoSetList"),
	  function(object, ...){
		  nms <- ls(assayData(object))
		  if("baf" %in% nms){
			  hmmBeadStudioSetList(object, ...)
		  } else hmmOligoSetList(object, ...)
	  })





