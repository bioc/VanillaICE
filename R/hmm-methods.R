setMethod("hmm", signature(object="oligoSnpSet"),
	  function(object, ...){
		  if("baf" %in% assayDataElementNames(object)){
			  return(hmmBafLrrSet2(object, ...))
		  } else{
			  ## baf not in assay data
			  ## - use genotypes
			  return(hmmOligoSnpSet2(object, ...))
		  }
	  })

setMethod("hmm", signature(object="BeadStudioSet"),
	  function(object, ...){
		  hmmBafLrrSet2(object, ...)
	  })

setMethod("hmm", signature(object="SnpSet2"),
	  function(object, ICE=FALSE, ...){
		  if(!ICE){
			  hmmSnpSet2(object, ...)
		  } else hmmSnpSetIce(object, ...)
	  })

setMethod("hmm", signature(object="BeadStudioSetList"),
	  function(object, ...){
		  ##hmmBeadStudioSetList(object, ...)
		  hmmBafLrrSetList2(object, ...)
	  })

setMethod("hmm", signature(object="BafLrrSetList"),
	  function(object, ...){
		  hmmBafLrrSetList2(object, ...)
	  })

setMethod("hmm", signature(object="oligoSetList"),
	  function(object, ...){
		  nms <- ls(assayData(object))
		  if("baf" %in% nms){
			  hmmBafLrrSetList2(object, ...)
		  } else {
			  stop("assay data element 'baf' required")
		  }##hmmOligoSetList(object, ...)
	  })





