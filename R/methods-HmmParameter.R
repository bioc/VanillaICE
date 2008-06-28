setMethod("hmm", c("HmmOptions", "HmmParameter"),
          function(object, params){
		  if(!validObject(object)){
			  stop("the object of class HmmOptions must be valid. See validObject()")
		  }
		  if(!validObject(params)){
			  stop("the object of class HmmParameter must be valid. See validObject()")
		  }
		  ##harmonize the two objects (can't do this without featureData in params)
		  ##i <- match(featureNames(object), featureNames(params))
		  ##j <- match(sampleNames(object), sampleNames(hmmPredict))
		  snpset <- snpset(object)
		  if("copyNumber" %in% assayDataElementNames(snpset) | "ratio" %in% assayDataElementNames(snpset)){		  
			  if(min(copyNumber(snpset), na.rm=TRUE) > 0){
				  print("Transforming copy number to log2 scale.")
				  copyNumber(snpset) <- log2(copyNumber(snpset))
			  } else{
				  print("Negative values in the copy number.  Assume that copy number has been suitably transformed and is approximately Gaussian")
			  }
		  }

		  ##Useful for trios, but otherwise should be ignored
		  if(dim(emission(params))[2] != length(sampleNames(snpset))){
			  sn <- dimnames(emission(params))[[2]]
			  if(is.null(sn)){
				  warning("dimnames not provided for emission probabilities")
				  sn <- "sample1"
			  }
		  } else{
			  sn <- sampleNames(snpset)
		  }
		  arm <- paste(chromosome(snpset), featureData(snpset)$arm, sep="")
		  beta <- emission(params)
		  predictions <- matrix(NA, nrow=nrow(snpset), ncol=(dim(beta)[2]))
		  for(i in 1:(dim(beta)[2])){
			  if(i == 1){
				  print(paste("Fitting HMM to sample", i))
			  } else {cat(i, "\n")}
			  predictions[, i] <- viterbi(beta[, i, ],
						      pi=pi(params),
						      tau=genomicDistance(params),
						      states=states(object),
						      arm=arm,
						      tau.scale=transitionScale(params))
		  }
		  rownames(predictions) <- featureNames(snpset)
		  colnames(predictions) <- sn
		  i <- match(sn, sampleNames(snpset))
		  hmmOut <- new("HmmPredict",
				states=states(object),
				predictions=predictions,
				featureData=featureData(snpset),
				experimentData=experimentData(snpset),
				phenoData=phenoData(snpset[, i]),
				annotation=annotation(snpset))
		  breakpoints(hmmOut) <- calculateBreakpoints(hmmOut)
		  return(hmmOut)
	  })

setMethod("emission", "HmmParameter", function(object) object@emission)
setReplaceMethod("emission", c("HmmParameter", "array"),
		 function(object, value){
			 object@emission <- value
			 object
		 })
setMethod("tau", "HmmParameter", function(object) object@genomicDistance)
setMethod("genomicDistance", "HmmParameter", function(object) object@genomicDistance)
setReplaceMethod("genomicDistance", c("HmmParameter", "numeric"),
		 function(object, value) {
			 object@genomicDistance <- value
			 object
			 })
setMethod("pi", "HmmParameter", function(object) object@initialStateProbability)
setMethod("transitionScale", "HmmParameter", function(object) object@transitionScale)
##setMethod("tau.scale", "HmmParameter", function(object) object@transitionScale)
setReplaceMethod("transitionScale", c("HmmParameter", "matrix"),
		 function(object, value){
			 object@transitionScale <- value
			 object
			 })

##setMethod("hmmOptions", "HmmParameter", function(object) object@hmmOptions)
##setReplaceMethod("hmmOptions", c("HmmParameter", "HmmOptions"),
##                 function(object, value){
##                   object@hmmOptions <- value
##                   object
##                 })
##setMethod("featureData", "HmmParameter", function(object) object@featureData)
##setMethod("featureNames", "HmmParameter", function(object) featureNames(featureData(object)))
##setMethod("nrow", "HmmParameter", function(x) length(featureNames(x)))
##setMethod("fData", "HmmParameter", function(object) pData(featureData(object)))
##setMethod("sampleNames", "HmmParameter", function(object) colnames(Beta(object)))
##setMethod("chromosome", "HmmParameter", function(object) featureData(object)$chromosome)
##setMethod("arm", "HmmParameter", function(object) featureData(object)$arm)
##setMethod("Beta", "HmmParameter", function(object) object@beta)
##setMethod("show", "HmmParameter", function(object) str(object))
##setMethod("SnpClass", "HmmParameter", function(object) SnpClass(hmmOptions(object)))
setMethod("states", "HmmParameter", function(object) object@states)
setReplaceMethod("states", c("HmmParameter", "character"),
                 function(object, value){
                   hmmOptions(object)@states <- value
                   object
                 })
##setMethod("pi", "HmmParameter", function(object) object@pi)
##setMethod("tau", "HmmParameter", function(object) object@tau)
##setMethod("gt.ICE", "HmmParameter", function(object) gt.ICE(hmmOptions(object)))
##setMethod("cn.ICE", "HmmParameter", function(object) cn.ICE(hmmOptions(object)))
##setMethod("cn.location", "HmmParameter", function(object) cn.location(hmmOptions(object)))
##setMethod("cn.robustSE", "HmmParameter", function(object) cn.robustSE(hmmOptions(object)))
##setMethod("tau.scale", "HmmParameter", function(object) object@tau.scale)
##setMethod("[", "HmmParameter",
##          function(x, i, j, ..., drop = FALSE){
##            if (missing(drop)) drop <- FALSE
##            if (missing(i) && missing(j)) {
##              if (length(list(...))!=0)
##                stop("specify genes or samples to subset; use '",
##                     substitute(x), "$", names(list(...))[[1]],
##                     "' to access phenoData variables")
##              return(x)
##            }
##            ##number of rows (SNPs)
##            R <- length(featureNames(x))
##            ##number of states 
##            S <- length(pi(x))
##            ##number of columns (samples)
##            C <- ncol(Beta(x))
##            if(!missing(i) & !missing(j)){
##              x@beta <- Beta(x)[i, , j, drop=FALSE]
##              x@featureData <- featureData(x)[i, ]
##	      idx <- sort(unique(c(i-1, i, i+1)))
##	      idx <- idx[idx != 0]
##	      idx <- idx[idx != R]
##	      x@tau <- tau(x)[idx]
##            }
##            if(!missing(i) & missing(j)){
##              x@beta <- Beta(x)[i, , , drop=FALSE]
##	      idx <- sort(unique(c(i-1, i, i+1)))
##	      idx <- idx[idx != 0]
##	      idx <- idx[idx != R]
##	      x@tau <- tau(x)[idx]
##              x@featureData <- featureData(x)[i, ]
##            }
##            x
##          })

setMethod("[", "HmmParameter",
          function(x, i, j, ..., drop = FALSE){
            if (missing(drop)) drop <- FALSE
            if (missing(i) && missing(j)) {
              if (length(list(...))!=0)
                stop("specify genes or samples to subset; use '",
                     substitute(x), "$", names(list(...))[[1]],
                     "' to access phenoData variables")
              return(x)
            }
            ##number of rows (SNPs)
            R <- dim(emission(x))[1]
            ##number of columns (samples)
            C <- dim(emission(x))[2]	    
            ##number of states 
            S <- length(states(x))
            if(!missing(i) & !missing(j)){
		    emission(x) <- emission(x)[i, j, , drop=FALSE]
		    idx <- sort(unique(c(i-1, i, i+1)))
		    idx <- idx[idx != 0]
		    idx <- idx[idx != R]
		    genomicDistance(x) <- genomicDistance(x)[idx]
            }
            if(!missing(i) & missing(j)){
		    emission(x) <- emission(x)[i, , drop=FALSE]
		    idx <- sort(unique(c(i-1, i, i+1)))
		    idx <- idx[idx != 0]
		    idx <- idx[idx != R]
		    genomicDistance(x) <- genomicDistance(x)[idx]
            }
	    if(missing(i) & !missing(j)){
		    emission(x) <- emission(x)[, j, drop=FALSE]
	    }
            x
          })

setMethod("show", "HmmParameter", function(object){
	str(object)
})

          
                   
