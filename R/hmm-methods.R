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

setMethod("hmm", signature(object="CNSet"),
	  function(object, ...){
		  hmmCNSet(object, ...)
	  })

setMethod("hmm", signature(object="SnpSet"),
	  function(object, ...){
		  hmmSnpSet(object, ...)
	  })

