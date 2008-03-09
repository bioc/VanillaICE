setMethod("hmm", "HmmParameter",
          function(object, snpset, sn, verbose=FALSE, ...){
		  if("copyNumber" %in% assayDataElementNames(snpset)){		  
			  if(min(copyNumber(snpset), na.rm=TRUE) > 0){
				  print("Transforming copy number to log2 scale.")
				  copyNumber(snpset) <- log2(copyNumber(snpset))
			  } else{
				  print("Negative values in the copy number.  Assume that copy number has been suitably transformed and is approximately Gaussian")
			  }
		  }
		  
		  arm <- paste(chromosome(object), arm(object), sep="")
		  beta <- Beta(object)
		  predictions <- matrix(NA, nrow=nrow(snpset), ncol=(dim(beta)[3]))
		  for(i in 1:(dim(beta)[3])){
			  if(verbose){
				  if(i == 1){
					  print(paste("Fitting HMM to sample", i))
				  } else {cat(i, "\n")}
			  }
			  predictions[, i] <- viterbi(beta[, , i],
						      pi=pi(object),
						      tau=tau(object),
						      states=states(object),
						      arm=arm,
						      tau.scale=tau.scale(object))
		  }
		  rownames(predictions) <- featureNames(object)
		  if(missing(sn)) sn <- sampleNames(snpset)
		  colnames(predictions) <- sn

		  hmmOut <- new("HmmPredict",
				states=states(object),
				predictions=predictions,
				SnpClass=SnpClass(object),
				featureData=featureData(object),
				experimentData=experimentData(snpset),
				phenoData=phenoData(snpset),
				annotation=annotation(snpset))
		  breakpoints(hmmOut) <- calculateBreakpoints(hmmOut)
		  return(hmmOut)
          })



setMethod("hmmOptions", "HmmParameter", function(object) object@hmmOptions)
setReplaceMethod("hmmOptions", c("HmmParameter", "HmmOptions"),
                 function(object, value){
                   object@hmmOptions <- value
                   object
                 })
setMethod("featureData", "HmmParameter", function(object) object@featureData)
setMethod("featureNames", "HmmParameter", function(object) featureNames(featureData(object)))
setMethod("nrow", "HmmParameter", function(x) length(featureNames(x)))
setMethod("fData", "HmmParameter", function(object) pData(featureData(object)))
setMethod("sampleNames", "HmmParameter", function(object) colnames(Beta(object)))
setMethod("chromosome", "HmmParameter", function(object) featureData(object)$chromosome)
setMethod("arm", "HmmParameter", function(object) featureData(object)$arm)
setMethod("Beta", "HmmParameter", function(object) object@beta)
setMethod("show", "HmmParameter", function(object) str(object))
setMethod("SnpClass", "HmmParameter", function(object) SnpClass(hmmOptions(object)))
setMethod("states", "HmmParameter", function(object) states(hmmOptions(object)))
setReplaceMethod("states", c("HmmParameter", "character"),
                 function(object, value){
                   hmmOptions(object)@states <- value
                   object
                 })
setMethod("pi", "HmmParameter", function(object) object@pi)
setMethod("tau", "HmmParameter", function(object) object@tau)
setMethod("gt.ICE", "HmmParameter", function(object) gt.ICE(hmmOptions(object)))
setMethod("cn.ICE", "HmmParameter", function(object) cn.ICE(hmmOptions(object)))
setMethod("cn.location", "HmmParameter", function(object) cn.location(hmmOptions(object)))
setMethod("cn.robustSE", "HmmParameter", function(object) cn.robustSE(hmmOptions(object)))
setMethod("tau.scale", "HmmParameter", function(object) object@tau.scale)

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
            R <- length(featureNames(x))
            ##number of states 
            S <- length(pi(x))
            ##number of columns (samples)
            C <- ncol(Beta(x))
            if(!missing(i) & !missing(j)){
              x@beta <- Beta(x)[i, , j, drop=FALSE]
              x@featureData <- featureData(x)[i, ]
	      idx <- sort(unique(c(i-1, i, i+1)))
	      idx <- idx[idx != 0]
	      idx <- idx[idx != R]
	      x@tau <- tau(x)[idx]
            }
            if(!missing(i) & missing(j)){
              x@beta <- Beta(x)[i, , , drop=FALSE]
	      idx <- sort(unique(c(i-1, i, i+1)))
	      idx <- idx[idx != 0]
	      idx <- idx[idx != R]
	      x@tau <- tau(x)[idx]
              x@featureData <- featureData(x)[i, ]
            }
            x
          })

          
                   
