calls.emission <- function(object){
	if(!("calls" %in% assayDataElementNames(object@snpset))){
		warning("calls not an element of the assayData for object@snpset.")
		return()
	}
	if(!validObject(object)){
		stop("HmmOptions object is not valid")
	}
	snpset <- object@snpset
	states <- object@states
	S <- length(states)	
	gte <- array(calls(snpset), c(nrow(snpset), ncol(snpset), S))
	dimnames(gte) <- list(featureNames(snpset), sampleNames(snpset), states)
	if(!object@calls.ICE) {
		##for algorithms that produce bi-allelic calls, only 3
		##calls are possible.  The emission probability for
		##the call is the emission probability for a
		##homozygous call versus the emission probability for
		##a heterozygous call

		i <- which(as.vector(gte) == 1 | as.vector(gte) == 3)
		j <- which(as.vector(gte) == 2)
		k <- which(as.vector(gte) == 4 | is.na(as.vector(gte)))

		##prob. for each state.
		cn.location <- aperm(array(object@copyNumber.location, dim=c(S, ncol(snpset), nrow(snpset))))
		probHom <-  aperm(array(object@probHomCall, dim=c(S, ncol(snpset), nrow(snpset))))
		
		emission.gt <- vector("numeric", length=length(as.vector(gte)))
		emission.gt[i] <- as.vector(probHom)[i]
		emission.gt[j] <- (1-as.vector(probHom))[j]
		emission.gt[k] <- NA
		emission.gt[c(i, j)] <- log(emission.gt[c(i, j)])
		emission.gt <- array(emission.gt, dim=dim(gte))
	}  else {
		pkgs <- strsplit(annotation(object@snpset), ",")[[1]]
		if(length(pkgs) > 1){
			object@snpset <- .splitAnnotation(snpset)
		} else{
			featureData(object@snpset)$platform <- NA
			featureData(object@snpset)$platform <- annotation(snpset)
		}
		emission.gt <- array(NA, dim=c(nrow(object@snpset), ncol(object@snpset), 2))
		dimnames(emission.gt) <- list(featureNames(object@snpset), sampleNames(object@snpset), c("LOH", "N"))
		for(i in 1:length(pkgs)){
			j <- which(fData(object@snpset)$platform == pkgs[i])
			annotation(object@snpset) <- pkgs[1]
			emission.gt[j, , ] <- .getCallEmission(object[j, ])
		}
		emission.gt <- log(emission.gt)
	}
	return(emission.gt)##log scale
}

copyNumber.emission <- function(object){
	if(!any(c("copyNumber", "ratio") %in% assayDataElementNames(object@snpset))){
		warning("'copyNumber' or 'ratio' must be elements of the assayData for object@snpset.")
		return()
	}
	snpset <- snpset(object)
	cne <- copyNumber(snpset)
	if(all(cne > 0, na.rm=TRUE)){
		print("Calculating emission probabilities on the log(copy number)")
		i <- which(!is.na(as.vector(cne)))

		##Take log of hidden state
		if(any(object@copyNumber.location < 0)){
			stop("copyNumber.location must be positive")
		}
		tmp <- rep(NA, length(as.vector(cne)))
		tmp[i] <- log2(cne[i])
		cne <- matrix(tmp, nrow=nrow(snpset), ncol=ncol(snpset))		
		cn.location <- log2(object@copyNumber.location)
		takelog <- TRUE
	} else{
		print("Negative values in copy number --  assume that copy number has been transformed to log scale")
		cn.location <- object@copyNumber.location
		takelog <- FALSE
	}
	fn <- names(cne)
	states <- object@states
	S <- length(states)

	cne <- array(cne, dim=c(nrow(cne), ncol(cne), S))
	dimnames(cne) <- list(featureNames(snpset), sampleNames(snpset), states)
	##assume true copy number mean is the same for all samples
	##Easiest to keep everything the same dimension, even if not needed at this point
	cn.location <- aperm(array(cn.location, dim=c(S, ncol(snpset), nrow(snpset))))
	
	if(!object@copyNumber.ICE){
		##Use robust estimate of standard error-- sample-specific
		i <- which(chromosome(snpset) %in% as.character(1:22))
		snpset.autosomes <- snpset[i, ]
		f <- function(x){
			diff(quantile(x, probs=c(0.16, 0.84), na.rm=TRUE))/2 ##+ quantile(x, probs=0.16)
		}
		cne.auto <- copyNumber(snpset.autosomes)
		if(takelog) {
			cne.auto <- log2(cne.auto)
		}

		if(is.null(object@copyNumber.scale)){
			cn.robustSE <- apply(cne.auto, 2, f)			
			copyNumber.scale <- array(matrix(cn.robustSE,
							 nrow=nrow(snpset),
							 ncol=ncol(snpset),
							 byrow=TRUE), dim=dim(cne))
		} else{
			copyNumber.scale <- aperm(array(object@copyNumber.scale, dim=c(S, ncol(snpset), nrow(snpset))))
			if(takelog) print("User-supplied copyNumber.scale should be a standard deviation of the log2 CN")
		} 
		dimnames(copyNumber.scale) <- dimnames(cne)		
		i <- which(!is.na(as.vector(cne)))
	} else{
		##Use SNP-specific standard errors
		i <- which(!(is.na(as.vector(cnConfidence(snpset)))) & !is.na(as.vector(cne)))
		if(all(cnConfidence(snpset) > 0, na.rm=TRUE)){
			print("Using 1/cnConfidence(object) as standard errors for the copy number")
			copyNumber.scale <- array(1/cnConfidence(snpset), dim=dim(cne))
		} else{
			stop("confidence scores in slot cnConfidence must be positive")
		}
		
	}
	k <- which(!is.na(as.vector(cne)))
	if(!identical(dim(cne), dim(cn.location))) stop("dimensions must be the same")
	emission.cn <- rep(NA, length(as.vector(cn.location)))
	emission.cn[k] <- log(dnorm(as.vector(cne)[k], as.vector(cn.location)[k], as.vector(copyNumber.scale)[k]))
	emission.cn <- array(emission.cn, dim=dim(cne))
	dimnames(emission.cn) <- dimnames(cne)
	return(emission.cn) ##log scale
}
				    

##setMethod("gtEmission", c("numeric", "HmmOptions"),
##          function(x, options, confidence, featureNames, ...){
##            gtEmission(as.integer(x), options, confidence, featureNames, ...)
##          })

##setMethod("gtEmission", c("integer", "HmmOptions"),
##          function(x, options, confidence, featureNames, ...){
##            gte <- x; rm(x)
##            fn <- featureNames
##            if(options@SnpClass == "SnpCopyNumberSet") return()
##            states <- options@states
##            S <- length(states)
##            
##            ##gte: genotype estimates
##            ##gte <- as.vector(gte)
##            gte[gte == 3] <- 1
##            gte[gte == 4 | is.na(gte)] <- 3 ##Make missing genotypes have value 3
##
##            if(!options@gt.ICE) {
##              gte.state <- options@gte.state
##              gte.state <- t(cbind(gte.state, 1-gte.state, NA))
##              gt.emission <- gte.state[gte, ]
##            }  else {
##              gte <- as.integer(gte)
##              names(gte) <- fn
##              gt.emission <- .gtEmission.ICE(x=gte,
##                                             confidence=confidence,
##                                             options=options)
##            }
##            log(gt.emission)
##          })



##If there is more than one sample in object, it uses only the first.
##setMethod("cnEmission", c("numeric", "HmmOptions"),
##          function(x, options, confidence, robustSE, ...){
##		  if(options@SnpClass == "SnpCallSet") return()            
##		  cne <- x; rm(x)
##		  fn <- names(cne)
##		  states <- options@states
##		  S <- length(states)
##		  cn.location <- options@cn.location
##		  cne <- matrix(cne, nrow=length(cne), ncol=S, byrow=FALSE)
##            
##		  ##assume true copy number mean is the same for all samples
##		  cn.location <- matrix(cn.location, nrow=nrow(cne), ncol=S, byrow=TRUE)
##
##		  if(!options@cn.ICE){
##			  ##Use robust estimate of standard error-- sample-specific
##			  cn.SE <- matrix(robustSE, nrow=nrow(cn.location), ncol=ncol(cn.location), byrow=TRUE)
##		  } else{
##			  ##Use SNP-specific standard errors
##			  cn.SE <- 1/confidence
##			  cn.SE <- matrix(cn.SE, nrow=nrow(cne), ncol=ncol(cne), byrow=FALSE)
##			  colnames(cn.SE) <- states
##			  if(any(is.na(cn.SE))) stop("NA's in confidence scores.  Must  exclude these SNPs or plug in values for the confidence scores")
##		  }
##		  emission.cn <- dnorm(cne, cn.location, cn.SE)
##		  ##rownames(emission.cn) <- names(cne)
##		  ##colnames(emission.cn) <- states
##		  ##length should be R*S (R = number of SNPs, S = number of states)
##		  log(emission.cn)
##          })
