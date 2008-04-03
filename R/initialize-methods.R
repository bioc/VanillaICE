

setMethod("initialize", "HmmParameter",
	  function(.Object,
		   states=character(),
		   initialStateProbability=0.99, ##probability of baseline state					   
		   emission=array(),
		   genomicDistance=numeric(),
		   transitionScale=matrix()){
		  .Object@states <- states
		  S <- length(states)
		  if(length(initialStateProbability) == 1){
			  pi <- initialStateProbability
			  pi.other <- (1-pi)/(length(states)-1)
			  pi <- rep(pi, length(states))
			  pi[states != "N"] <- pi.other
			  initialStateProbability <- pi
		  }
		  .Object@initialStateProbability <- initialStateProbability
		  .Object@emission <- emission
		  .Object@genomicDistance <- genomicDistance
		  .Object@transitionScale <- transitionScale
		  .Object
	  })
setValidity("HmmParameter",
	    function(object){
		    valid <- TRUE
		    msg <- NULL
		    S <- length(object@states)
		    if(S != length(pi(object))){
			    msg <- c(msg, "length of states not equal to length of initial state probabilities")
			    valid <- FALSE
		    }
		    if(S != ncol(object@transitionScale) | S != nrow(object@transitionScale)){
			    msg <- c(msg, "transitionScale must have the same number of rows and columns as there are hidden states")
			    valid <- FALSE
		    }
		    if(S != dim(emission(object))[3]){
			    msg <- c(msg, "the third dimension for the array of emission probabilities must be equal to the number of hidden states")
			    valid <- FALSE
		    }
		    if(length(object@genomicDistance) != (dim(emission(object))[1]-1)){
			    msg <- c(msg, "genomicDistance should be a numeric vector of length = number of SNPs - 1")
			    valid <- FALSE
		    }
		    if(!valid) print(msg)
		    valid
	    })
		       
			    


setMethod("initialize", "HmmOptions",
	  function(.Object,
		   snpset,
		   states="",
		   copyNumber.location=numeric(),
		   copyNumber.scale=NULL,
		   copyNumber.ICE=FALSE,
		   calls.ICE=FALSE,
		   probHomCall=vector("numeric", length(states)),
		   probWrongCall=NULL){
		  if(missing(snpset)) stop("snpset must be specified")
		  if(any(is.na(chromosome(snpset)))){
			  snpset <- snpset[!is.na(chromosome(snpset)), ]
		  }
		  if(!all(chromosome(snpset) %in% c(1:22, "X"))){
			  print("Excluding SNPs that are not on chromosomes 1-22 or X")
			  snpset <- snpset[chromosome(snpset) %in% c(1:22, "X"), ]
		  }
		  chrom <- chromosome2numeric(chromosome(snpset))
		  if(class(snpset) == "RatioSnpSet"){
			  featureAnn <- which(!c("chromosome", "position", "arm") %in% fvarLabels(snpset))
			  if(length(featureAnn) > 0){
				  stop(paste(c("chromosome", "position", "arm")[featureAnn], "must be added to the featureData"))
			  }
		  }
		  if(!("chromosome" %in% fvarLabels(snpset))){
			  featureData(snpset)$chromosome <- chromosome(snpset)
		  }
		  if(!("position" %in% fvarLabels(snpset))){
			  featureData(snpset)$position <- position(snpset)
		  }
		  if(!("arm" %in% fvarLabels(snpset))){
			  featureData(snpset)$arm <- getChromosomeArm(snpset)
			  if(any(is.na(featureData(snpset)$arm))) warning("NAs for chromosomal arm")
		  }
		  .Object@snpset <- snpset[order(chrom, position(snpset)), ]
		  .Object@states <- states
		  .Object@copyNumber.location <- copyNumber.location
		  .Object@copyNumber.scale <- copyNumber.scale
		  .Object@copyNumber.ICE <- copyNumber.ICE
		  .Object@calls.ICE <- calls.ICE
		  .Object@probHomCall <- probHomCall
		  .Object@probWrongCall <- probWrongCall
		  .Object
	  })

setValidity("HmmOptions", function(object){
	##must specify the probability of a wrong call. order must be
	##for calls=1 (AA), calls=2 (AB), calls=3(BB)
	msg <- NULL
	valid <- TRUE
	if(!extends(class(object@snpset), "SnpLevelSet")){
		msg <- c(msg, "snpset must extend SnpLevelSet")
		valid <- FALSE
	}
	if(!(all(c("position", "arm", "chromosome") %in% fvarLabels(object@snpset)))){
		msg <- c(msg, "snpset must have chromosome, physical position, and (chromosomal) arm stored in the featureData.  See addFeatureData -- may require pd.mapping packages")
		valid <- FALSE
	}
	if(any(is.na(chromosome(object@snpset)))){
		msg <- c(msg, "snpset can not have SNPs with NA for chromosome")
		valid <- FALSE
	}
	if(!all(chromosome(object@snpset) %in% c(1:22, "X"))){
		msg <- c(msg, "Object contains chromosomes that are not currently supported for the HMM.  Supported chromosomes are 1-22 and X.")
		valid <- object
	}
	chrom <- chromosome2numeric(chromosome(object@snpset))
	tmp <- object@snpset[order(chrom, position(object@snpset)), ]
	if(!identical(featureNames(tmp), featureNames(object@snpset))){
		msg <- c(msg, "SNPs must be ordered by chromosome and physical position")
		valid <- FALSE
	}
	if(object@calls.ICE) {
		if(is.null(object@probWrongCall)){
			msg <- c(msg, "probability of a wrong genotype call must be specified when calls.ICE is TRUE. probWrongCall must be a numeric vector of length 3. Value in interval 0, 1.  Order is for calls=1 (AA), calls=2 (AB), calls=3(BB)")
			valid <- FALSE
		} else{
			if(length(object@probWrongCall) != 3){
				msg <- c(msg, "probWrongCall must be a numeric vector of length 3 with values in the interval 0, 1.  Order is for calls=1 (AA), calls=2 (AB), calls=3(BB)")
				valid <- FALSE
			}
		}
	}
	if(any(c("copyNumber", "ratio") %in% assayDataElementNames(object@snpset))){
		if(length(object@states) != length(object@copyNumber.location)){
			msg <- c(msg, "copyNumber.location (a numeric vector) must be specified for each hidden state (and in the same order).")
			valid <- FALSE
		}
		if(!is.null(object@copyNumber.scale)){
			if(!(all(dim(object@copyNumber.scale) == dim(copyNumber(object@snpset))))){
				valid <- FALSE
				msg <- "scale parameter for copy number must be the same dimension as the copy number matrix"
			}
		}
	}
	if("calls" %in% assayDataElementNames(object@snpset)){
		if(length(object@states) != length(object@probHomCall)){
			msg <- c(msg, "the probability of observing a homozygous call given the hidden state ('probHomCall') must be a numeric vector) must be specified for each hidden state (and in the same order).")
			valid <- FALSE
		}
		if(max(object@probHomCall) > 1 | min(object@probHomCall) < 0){
			msg <- c(msg, "values of probHomCall must be in the interval 0, 1")
			valid <- FALSE
		}		
	}
	if(!valid) print(msg)
	valid
})


##One big matrix that could potentially span several chromosomes.
##featureData should be looked up in RSQLite database
##setMethod("initialize", "HmmParameter",
##          function(.Object,
##                   snpset,
##                   states=character(),
##                   tau=numeric(),
##                   tau.scale=matrix(),                   
##                   pi=numeric(),
##                   notes=character(),
##                   beta=array(),
##                   SCALE=2,##by default, twice as likely to return to normal state
##		   ...){
##		  if(any(is.na(chromosome(snpset)))){
##			  warning("missing values for chromosome(snpset) -- these SNPs are removed")
##			  snpset <- snpset[!(is.na(chromosome(snpset))), ]			  
##		  }
##
##		  chrom <- chromosome2numeric(chromosome(snpset))
##		  tmp <- snpset[order(chrom, position(snpset)), ]
##		  if(!identical(featureNames(tmp), featureNames(snpset))){
##			  warning("snpset object ordered by chromosome and physical position.  Results will be returned in this order")
##		  }
##		  snpset <- tmp; rm(tmp)
##		  gc()
##
##		  if("copyNumber" %in% assayDataElementNames(snpset)){
##			  if(min(copyNumber(snpset), na.rm=TRUE) > 0){
##				  print("Transforming copy number to log2 scale.")
##				  copyNumber(snpset) <- log2(copyNumber(snpset))
##				  log2cn.location <- TRUE
##			  } else{
##				  print("Negative values in the copy number.  Assume that copy number has been suitably transformed and is approximately Gaussian")
##				  log2cn.location <- FALSE
##			  }
##		  }
##		  hmmOptions <- new("HmmOptions", snpset=snpset,
##				    states=states, beta=beta, ...)
##
##		  if("copyNumber" %in% assayDataElementNames(snpset)){		  
##			  if(log2cn.location){
##				  hmmOptions@cn.location <- log2(hmmOptions@cn.location)
##			  }
##		  }
##		  S <- length(hmmOptions@states)
##		  featureNames <- featureNames(snpset)
##		  featureData <- addFeatureData(snpset)
##
##                  ###########################################################################            
##                  ##Initial State Probabilities
##                  ###########################################################################
##		  states <- hmmOptions@states
##		  if(length(pi) == 0){
##			  tmp <- rep(0, length(states))
##			  names(tmp) <- states
##			  tmp["N"] <- 0.999
##			  tmp[states != "N"] <- (1-0.999)/(length(states) - 1)
##			  pi <- log(tmp)
##		  }
##
##                  ###########################################################################            
##		  ##Emission Probabilities
##                  ###########################################################################
##		  ##Rows = number of SNPs * number of states, Columns = number of samples
##		  ##array: dims = number of SNPs, number states, number samples
##		  if(length(beta) == 1){
##			  if(hmmOptions@SnpClass == "SnpCallSet" | hmmOptions@SnpClass == "oligoSnpSet"){
##				  x <- data.frame(calls(snpset), row.names=featureNames(snpset))
##				  confidence <- data.frame(callsConfidence(snpset), row.names=featureNames(snpset))
##				  print("Calculating emission probabilities for genotype calls.... ")
##				  beta.gt <- array(NA, c(nrow(snpset), length(states), ncol(snpset)))
##				  for(i in 1:ncol(snpset)){
##					  beta.gt[, , i] <- gtEmission(x=x[, i],
##								       confidence=confidence[, i],
##								       options=hmmOptions,
##								       featureNames=featureNames(snpset))
##				  }
##			  } else { beta.gt <- 0}
##			  if(hmmOptions@SnpClass == "SnpCopyNumberSet" | hmmOptions@SnpClass == "oligoSnpSet"){
##				  x <- data.frame(copyNumber(snpset), row.names=featureNames(snpset))
##				  confidence <- data.frame(cnConfidence(snpset), row.names=featureNames(snpset))
##				  cn.robustSE <- as.list(hmmOptions@cn.robustSE)
##				  print("Calculating emission probabilities for copy number estimates... ")
##				  beta.cn <- array(NA, c(nrow(snpset), length(states), ncol(snpset)))
##				  for(i in 1:ncol(snpset)){
##					  beta.cn[, , i] <- cnEmission(x=x[, i],
##								       confidence=confidence[, i],
##								       robustSE=cn.robustSE[[i]],
##								       options=hmmOptions,
##								       featureNames=featureNames(snpset))
##				  }
##			  } else {beta.cn <- 0}
##			  ##Missing values in the emission
##			  ##probabilities are typically caused by missing values for
##			  ##genotype calls.  Plugging in a constant such as zero for
##			  ##each of the hidden states should prevent the genotype
##			  ##call of the SNP from contributing to the likelihood,
##			  ##however the copy number estimate (if available) will
##			  ##still contribute
##			  beta.gt[is.na(beta.gt)] <- 0
##			  beta.cn[is.na(beta.cn)] <- 0
##			  beta <- beta.cn + beta.gt
##		  }
##
##                  ###########################################################################            
##		  ##Transition Probabilities
##                  ###########################################################################                        
##		  ##Scaling transition probabilities between states
##		  if(missing(tau.scale) & S > 2){
##			  ##scalar parameters for the transition probability
##			  ##matrix.  By default, the probability of transitioning
##			  ##from an altered state back to a normal state is twice
##			  ##as likely as the probability of transitioning between
##			  ##two altered states              
##			  S <- length(states)
##			  tau.scale <- matrix(1, S, S)
##			  rownames(tau.scale) <- colnames(tau.scale) <- states
##			  x <- seq(1, (S-1), by=0.01)
##			  y <- rep(NA, length(x))
##			  for(i in 1:length(x)){
##				  y[i] <- ((S-1)-x[i])/(S-2) 
##			  }
##			  z <- abs(x/y - SCALE)
##			  z <- z == min(z)
##			  x <- x[z]
##			  y <- y[z]
##			  tau.scale[upper.tri(tau.scale)] <- y
##			  tau.scale[lower.tri(tau.scale)] <- y
##			  tau.scale[, "N"] <- x
##			  ##by default, do not tau.scale the probabilities of leaving the normal state
##			  tau.scale["N", ] <- 1
##			  diag(tau.scale) <- 1
##		  }
##		  if(S == 2) tau.scale <- matrix(1, nrow=S, ncol=S)
##              
##		  ## We definine the transition probability matrix as
##		  ## a function of the distance between adjacent SNPs
##		  ## on the same chromosomal arm.  This requires
##		  ## breaking the object into a list by chromosome.
##		  ## We then order the SNPs by physical position along
##		  ## the chromosome, calculating the distance between
##		  ## adjacent SNPs. If there are S states, the
##		  ## transition probability matrix from SNP t to SNP
##		  ## t+1 has dimension S x S.  For each SNP t, we
##		  ## convert this matrix to a vector of length S^2.
##		  ## We then return a matrix of the transition
##		  ## probabilities with R-1 rows and S^2 columns. For
##		  ## SNPs at the beginning or the end of a chromosome,
##		  ## the transition probability matrix will be NA's.
##		  ## For SNPs at the beginning or end of a chromosomal
##		  ## arm, the transition probability matrix can be any
##		  ## number (it will not be used by the HMM).
##		  if(length(tau) == 0){
##			  ##              tau <- calculateTransitionProbability(object=snpset, options=hmmOptions)
##			  d <- calculateDistance(object=snpset)
##			  tau <- exp(-2 * d/(100*1e6))
##		  }
##		  if(any(is.na(tau))){
##			  warning("Missing values in transition probabilities")
##			  browser()
##		  }
##		  ##            names(tau) <- paste(featureNames(snpset)[1:(nrow(snpset)-1)], featureNames(snpset)[2:nrow(snpset)], sep="->")
##		  if("copyNumber" %in% assayDataElementNames(snpset)){		  
##			  if(log2cn.location){
##				  hmmOptions@cn.location <- 2^hmmOptions@cn.location
##			  }
##		  }
##		  .Object@hmmOptions <- hmmOptions
##		  .Object@tau <- tau
##		  .Object@tau.scale <- tau.scale
##		  .Object@pi <- pi
##		  .Object@beta <- beta
##		  .Object@featureData <- featureData
##		  .Object
##          })


setMethod("initialize", "HmmPredict",
          function(.Object,
                   predictions=new("matrix"),
		   SnpClass=character(),
		   states=character(),
		   breakpoints=data.frame(), ... ){
		  .Object <- callNextMethod(.Object,
					    predictions=predictions, ...)
##		  .Object@SnpClass <- SnpClass
		  .Object@states <- states
		  breakpoints(.Object) <- breakpoints
		  .Object
          })


##setMethod("initialize", "HmmPredict",
##          function(.Object,
##                   states=character(),
##                   predictions=matrix(),
##                   breakpoints=list(),
##                   SnpClass=character(),
##                   featureData=new("AnnotatedDataFrame")){
##            .Object@predictions <- predictions
##            .Object@breakpoints <- breakpoints
##            .Object@SnpClass <- SnpClass
##            .Object@states <- states
##            .Object@featureData <- featureData
##            .Object
##          })


setValidity("HmmPredict", function(object) {
  assayDataValidMembers(assayData(object), "predictions")
})




