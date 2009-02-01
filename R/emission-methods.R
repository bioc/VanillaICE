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
		##cn.location <- aperm(array(object@copyNumber.location, dim=c(S, ncol(snpset), nrow(snpset))))
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

copynumberEmission <- function(copynumber,
			       states,
			       mu,
			       sds,
			       takeLog,
			       verbose=TRUE){
	if(missing(sds)) stop("must supply sds estimates")		
	if(!all(sds > 0, na.rm=TRUE)) stop("sds not positive")
	cne <- copynumber
	if(missing(takeLog)) stop("must specify whether to take the log2 of the copy number matrix")
	if(takeLog){
		cne <- log2(cne)
		mu <- log2(mu)
	} 
	fn <- rownames(cne)
	S <- length(states)
	cne <- array(cne, dim=c(nrow(cne), ncol(cne), S))
	scale <- aperm(array(sds, dim=c(S, ncol(cne), nrow(cne))))
	##dimnames(cne) <- list(rownames(copynumber), colnames(copynumber), states)
	mu <- aperm(array(mu, dim=c(S, ncol(cne), nrow(cne))))
	##Use SNP-specific standard errors
	##i <- which(!(is.na(as.vector(1/sds))) & !is.na(as.vector(cne)))
	k <- which(!is.na(as.vector(cne)))
	if(!identical(dim(cne), dim(mu))) stop("dimensions must be the same")
	emission.cn <- rep(NA, length(as.vector(mu)))
	emission.cn[k] <- dnorm(as.vector(cne)[k], as.vector(mu)[k], as.vector(scale)[k])
	emission.cn <- array(emission.cn, dim=dim(cne))
	##dimnames(emission.cn) <- list(rownames(cne), colnames(cne), states)
	if(verbose) message("returning the emission probability on the log scale")
	return(log(emission.cn))
}

genotypeEmission <- function(genotypes, conf, states, probHomCall, probMissing, verbose=TRUE){
	##function to replace .getCallEmission
	if(!is.numeric(genotypes)) stop("genotypes must be integers (1=AA, 2=AB, 3=BB, 4=missing")
	emission <- array(genotypes, dim=c(nrow(genotypes), ncol(genotypes), length(states)))
	for(s in seq(along=states)){
		tmp <- genotypes
		tmp[tmp == 1 | tmp == 3] <- probHomCall[s]
		tmp[tmp == 2] <- 1-probHomCall[s]
		if(!missing(probMissing)) tmp[tmp == 4 | is.na(tmp)] <- probMissing[s]
		emission[, , s] <- tmp
	}
	dimnames(emission)[[3]] <- states
	logemit <- log(emission)
	return(logemit)
}

genotypeEmissionCrlmm <- function(genotypes, conf,
				  pHetCalledHom=0.001,
				  pHetCalledHet=0.995,
				  pHomInNormal=0.99,
				  pHomInLoh=0.999,
				  cdfName){
	##function to replace .getCallEmission
	if(missing(cdfName)) stop("must provide cdfName")
	if(cdfName != "GenomeWideSnp6") stop("currently this function only works for Affy 6.0")
	require(callsConfidence) || stop("callsConfidence package not available")

	if(cdfName == "GenomeWideSnp6"){
		data(affy6, package="callsConfidence", envir=.callsConfidencePkgEnv)
		hapmapP <- affy6
		hapmapP[, 2] <- 1-exp(-hapmapP[, 2]/1000)
		confidence <- 1-exp(-conf/1000)
	}
	i11 <- hapmapP[, 1] == 3  ##called homozygous truth homozygous
	i12 <- hapmapP[, 1] == 4  ##called homozygous truth heterozygous
	i21 <- hapmapP[, 1] == 1  ##called het truth hom
	i22 <- hapmapP[, 1] == 2  ##called het truth het
	f11 <- density(hapmapP[i11, 2], from=0, to=1, n=1e3)
	f12 <- density(hapmapP[i12, 2], from=0, to=1, n=1e3)
	f21 <- density(hapmapP[i21, 2], from=0, to=1, n=1e3)
	f22 <- density(hapmapP[i22, 2], from=0, to=1, n=1e3)

	##-------------------------------------------------------------------------
	##distribution of observed call probabilities when the call is homozygous
	##-------------------------------------------------------------------------

	##-------------------------------------------------------------------------
	##P(phat | LOH, gte)  
	##-------------------------------------------------------------------------	
	GT <- as.integer(genotypes)
	pTruthIsNormal <- pTruthIsLoh <- rep(NA, length(GT))	
	confidence <- as.numeric(confidence)
	confidence[confidence==0] <- 0.01 ##Otherwise, NA's result
	
	hom <- which(GT == 1 | GT == 3)
	observedPcalledHom <- cut(confidence[hom], breaks=f11$x, labels=FALSE)
	pTruthIsLoh[hom] <- f11$y[observedPcalledHom]

	het <- which(GT == 2)
	observedPcalledHet <- cut(confidence[het], breaks=f11$x, labels=FALSE)
	pTruthIsLoh[het] <- f21$y[observedPcalledHet]

	##-------------------------------------------------------------------------
	##Calculate P(phat | Normal, HOM)
	##-------------------------------------------------------------------------					
	chet1 <- f22$y[cut(confidence[het], breaks=f22$x, labels=FALSE)] 
	chet2 <- f21$y[cut(confidence[het], breaks=f21$x, labels=FALSE)]
	##term5[1]=P(true genotype is HET | genotype call is AB, state is normal)
	pTruthIsNormal[het] <- chet1*pHetCalledHet + chet2*(1-pHetCalledHet)
	##chom1=called homozygous truth heterozygous
	chom1 <- f12$y[cut(confidence[hom], breaks=f12$x, labels=FALSE)]
	##chom2=called homozygous truth homozygous
	chom2 <- f11$y[cut(confidence[hom], breaks=f11$x, labels=FALSE)]
	##chom4 <- 0.9999    ##P(HOM|CHOM)
	##probability that the true state is HOM when genotype call is homozygous
	##pHetCalledHom = P(true genotype is HET | calls is AA or BB, state is normal)
	pTruthIsNormal[hom] <- chom1*pHetCalledHom + chom2*(1-pHetCalledHom)

	fNormal <- fLoh <- rep(NA, length(GT))
	fNormal[hom] <- pHomInNormal * pTruthIsNormal[hom]
	fNormal[het] <- (1-pHomInNormal) * pTruthIsNormal[het]
	fLoh[hom] <- pHomInLoh * pTruthIsLoh[hom]
	fLoh[het] <- (1-pHomInLoh) * pTruthIsLoh[het]

	f <- array(NA, dim=c(nrow(genotypes), ncol(genotypes), 2))
	dimnames(f)[[3]] <- c("Normal", "LOH")
	f[, , "Normal"] <- matrix(fNormal, nrow(genotypes), ncol(genotypes))
	f[, , "LOH"] <- matrix(fLoh, nrow(genotypes), ncol(genotypes))
	f[f  == 0] <- min(f[f > 0], na.rm=TRUE)
	f <- log(f)
	return(f)	
}

.getCallEmission <- function(object){
	##load hapmap probabilities
	probHomCall <- object@probHomCall
	term5 <- object@term5
	object <- object@snpset
	.hapmapProbabilities <- function(annotationPackage){
		##require("callsConfidence") || stop("callsConfidence package not available")
		print("callsConfidence package must be available")
		get(switch(annotationPackage,
			   pd.mapping50k.hind240=data(hindPhat),
			   pd.mapping50k.xba240=data(xbaPhat),
			   pd.mapping250k.nsp=data(nspPhat), ##need to calculate phats for 500k
			   pd.mapping250k.sty=data(styPhat),
			   stop("annotation package not supported")))
	}            	
	hapmapP <- .hapmapProbabilities(annotation(object))
	
	##map observed call probablilities to a nonparametric estimate of
	##the density using the reference p distributions when the truth is
	##known
	gte <- array(calls(object), dim=c(nrow(object), ncol(object), 2))
	dimnames(gte) <- list(featureNames(object), sampleNames(object), c("LOH", "Normal"))
	
	confidence <- callsConfidence(object)
	i11 <- hapmapP[, 1] == 3  ##called homozygous truth homozygous
	i12 <- hapmapP[, 1] == 4  ##called homozygous truth heterozygous
	i21 <- hapmapP[, 1] == 1  ##called het truth hom
	i22 <- hapmapP[, 1] == 2  ##called het truth het
	f11 <- density(hapmapP[i11, 2], from=1/1e3, to=1, n=1e3)
	f12 <- density(hapmapP[i12, 2], from=1/1e3, to=1, n=1e3)
	f21 <- density(hapmapP[i21, 2], from=1/1e3, to=1, n=1e3)
	f22 <- density(hapmapP[i22, 2], from=1/1e3, to=1, n=1e3)
	##-------------------------------------------------------------------------
	##distribution of observed call probabilities when the call was homozygous
	##-------------------------------------------------------------------------	
	pTruthIsNormal <- pTruthIsLoh <- matrix(NA, nrow=nrow(object), ncol=ncol(object))
	for(i in 1:ncol(object)){
		##-------------------------------------------------------------------------
		##P(phat | LOH, gte)  
		##-------------------------------------------------------------------------				
		hom <- which(gte[, i, 1] == 1 | gte[, i, 1] == 3)
		observedPcalledHom <- cut(confidence[hom, i],
					  breaks=f11$x,
					  labels=FALSE) #called HOM, truth HOM w
		pTruthIsLoh[hom, i] <- f11$y[observedPcalledHom]
		
		het <- which(gte[, i, 1] == 2)
		observedPcalledHet <- cut(as.vector(confidence[het, i]),
					  breaks=f11$x,
					  labels=FALSE)
		pTruthIsLoh[het, i] <- f21$y[observedPcalledHet]		
		
		##-------------------------------------------------------------------------
		##Calculate P(phat | Normal, HOM)
		##-------------------------------------------------------------------------					
		chet1 <- f22$y[cut(confidence[het, i], breaks=f22$x, labels=FALSE)] 
		chet2 <- f21$y[cut(confidence[het, i], breaks=f21$x, labels=FALSE)]
		##term5[1]=P(true genotype is HET | calls is AB, state is normal)
		pTruthIsNormal[het, i] <- chet1*term5[1] + chet2*(1-term5[1])

		##chom1=called homozygous truth heterozygous
		chom1 <- f12$y[cut(confidence[hom, i], breaks=f12$x, labels=FALSE)]
		##chom2=called homozygous truth homozygous
		chom2 <- f11$y[cut(confidence[hom, i], breaks=f11$x, labels=FALSE)]
		##chom4 <- 0.9999    ##P(HOM|CHOM)
		##probability that the true state is HOM when genotype call is homozygous
		##term5[2]=P(true genotype is HET | calls is AA or BB, state is normal)
		pTruthIsNormal[hom, i] <- chom1*(term5[2]) + chom2*(1-term5[2])
	}

	stateProbability <- array(NA, dim=dim(gte))
	dimnames(stateProbability) <- dimnames(gte)
	stateProbability[, , 1] <- pTruthIsLoh
	stateProbability[, , 2] <- pTruthIsNormal
	B <- cbind(probHomCall, 1-probHomCall)
	colnames(B) <- c("AAorBB", "AB")
##	B <- cbind(gte.state, 1-gte.state, NA)
##	colnames(B) <- c("gte=HOM", "gte=HET", "gte=NA")
	rownames(B) <- c("LOH", "N")
##	rownames(B) <- options@states
	fobs <- array(NA, dim=c(nrow(object), length(sampleNames(object)), 2))
	dimnames(fobs) <- list(featureNames(object), sampleNames(object), c("LOH", "N"))
	for(i in 1:ncol(object)){
		hom <- which(gte[, i, 1] == "1" | gte[, i, 1] == "3")
		het <- which(gte[, i, 1] == "2")
		fobs[hom, i, ] <- matrix(B[, 1], nrow=length(hom), ncol=2, byrow=TRUE)
		fobs[het, i, ] <- matrix(B[, 2], nrow=length(het), ncol=2, byrow=TRUE)
	}
	emission.gt <- fobs * stateProbability
	dimnames(emission.gt) <- dimnames(gte)
	emission.gt
}
				    

