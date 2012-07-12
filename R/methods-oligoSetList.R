setMethod("[[", signature(x="oligoSetList"),
	  function(x, i, j, ..., exact=TRUE){
		  if(missing(i)) return(x)
		  if(length(i) == 1){
			  r <- copyNumber(x)[[i]]
			  b <- baf(x)[[i]]
			  gt <- calls(x)[[i]]
			  gtp <- snpCallProbability(x)[[i]]
			  colnames(r) <- colnames(b) <- colnames(gt) <- colnames(gtp) <- sampleNames(phenoData(x))
			  fd <- featureDataList(x)[[i]]
			  x <- new("oligoSnpSet",
				   call=gt,
				   callProbability=gtp,
				   copyNumber=r,
				   baf=b,
				   phenoData=phenoData(x),
				   featureData=fd,
				   genome=genomeBuild(x),
				   annotation=annotation(x))
		  } else {
			  stop("subscript out of bounds")
		  }
	  })


