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


setMethod("hmm", signature(object="BeadStudioSetList"),
	  function(object, ...){
		  message("Fitting HMM to each chromosome")
		  if(isPackageLoaded("ff")) pkgs <- c("ff", "VanillaICE") else pkgs <- "VanillaICE"
		  if(is.null(getCluster())) registerDoSEQ()
		  if(ncol(object[[1]]) < ocSamples()){
			  rdl <- foreach(obj=object, .packages=pkgs) %dopar% {
				  hmm(object=obj, ...)
			  }
			  ##rdl <- foreach(obj=object, .packages=pkgs, .combine="c") %dopar% nrow(obj)
		  } else {
			  sample.index <- splitIndicesByLength(seq_len(ncol(object[[1]])), ocSamples())
			  ## The outerloop is chromosome.  if the sample
			  ## size is very large, it probably makes
			  ## sense to have samples in the innerloop...
			  rdl <- foreach(obj=object) %:% foreach(j=sample.index, .packages=pkgs) %dopar% {
				  ## all samples, 1 chromosome is the inner-loop
				  hmm(object=obj[, j], ...)
			  }
			  ## outerloop is samples
##			  rdl <- foreach(j=sample.index) %:% foreach(obj=object, .packages=pkgs) %dopar% {
##				  ## the inner loop is chromosome
##				  ## (does all chromosomes for ocSamples())
##				  hmm(object=obj[, j], ...)
##			  }
			  rdl <- lapply(rdl, stackRangedData)
		  }
		  rd <- stackRangedData(rdl)
		  return(rd)
	  })

