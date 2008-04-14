##setMethod(".calls_ICE", c("integer", "HmmOptions"),
.splitAnnotation <- function(object, ...){  ##, P.CHOM.Normal, P.CHOM.LOH, SAMPLE=1){
	pkgs <- strsplit(annotation(object), ",")[[1]]	
	fn <- featureNames(object)
	if(length(pkgs) > 1){
		require(RSQLite) || stop("RSQLite package not available")
		sql <- "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
		object1 <- object2 <- object
		annotation(object1) <- pkgs[1]
		annotation(object2) <- pkgs[2]
		tmp1 <- dbGetQuery(db(object1), sql)
		tmp2 <- dbGetQuery(db(object2), sql)

		idx.pkg1 <- match(tmp1[["man_fsetid"]], fn)
		idx.pkg1 <- idx.pkg1[!is.na(idx.pkg1)]

		idx.pkg2 <- match(tmp2[["man_fsetid"]], fn)
		idx.pkg2 <- idx.pkg2[!is.na(idx.pkg2)]

		featureData(object)$platform <- rep(NA, nrow(object))
		featureData(object)$platform[idx.pkg1] <- pkgs[1]
		featureData(object)$platform[idx.pkg2] <- pkgs[2]
	}
	object
}
	
	
	
##		tmp1 <- .getCallEmission(object1)
##		tmp1 <- .getCallEmission(x=gte[enzyme == pkgs[1]],
##					 confidence=confidence[enzyme == pkgs[1]],
##					 hapmapP=hapmapP[[1]],
##					 options=options)
##		tmp2 <- .getCallEmission(x=gte[enzyme == pkgs[2]],
##					 confidence=confidence[enzyme == pkgs[2]],
##					 hapmapP=hapmapP[[2]],
##					 options=options)
##		tmp <- rbind(tmp1, tmp2)
##		idx <- match(fn, rownames(tmp))
##		gt.emission <- tmp[idx, ]
##		stopifnot(identical(rownames(gt.emission), fn))
##	} else{
##		gt.emission <- .getCallEmission(x=gte,
##						confidence=confidence,
##						hapmapP=hapmapP,
##						options=options)
##	} 
##	gt.emission
##}

##.getCallEmission <- function(x,  ##genotype estimates from oligo
##                             confidence, ##confidence scores from oligo
##                             hapmapP,  ##obtained from hapmap
##                             options){  ##HmmOptions object
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
##		observedPcalledHom <- cut(as.vector(confidence[gte == 1]),
##					       breaks=f11$x,
##					       labels=FALSE) #called HOM, truth HOM w
	##marginal distribution of observed call probabilities for heterozygous calls  
##	observedPcalledHet <- cut(as.vector(confidence[gte == 2]), breaks=f21$x,
##				  labels=FALSE) #called HET, truth HOM

	##this is approximately P(phat | HOM, gte=HET)
##	browser()
	##pTruthIsLoh <- rep(NA, length(gte))
##	pTruthIsLoh <- matrix(NA, nrow=nrow(object), ncol=ncol(object))
##	for(i in 1:ncol(object)){
##		hom <- which(gte[, i, 1] == "1" | gte[, i, 1] == "3")
##		pTruthIsLoh[hom, 1] <- f11$y[observedPcalledHom]
##	}


	##this is approximately P(phat | HOM, gte=HOM)  
##	pTruthIsLoh[gte == 2] <- f21$y[observedPcalledHet]

	##-------------------------------------------------------------------------
	##Calculate P(phat | Normal, CHET)
	##-------------------------------------------------------------------------	
	##This is approx. P(phat | HET, CHET) * P(HET | CHET) + P(phat | HOM, CHET) * P(HOM|CHET)
	##When the true state is normal, we need to estimate 4 things conditional on the call
	## P(p | CHET, normal):
	##chet1  P(p | HET, CHET)
	##chet2  P(p | HOM, CHET)
	##chet3  P(HET | CHET)
	##chet4  P(HOM | CHET)
##	pTruthIsNormal <- rep(NA, length(gte))
##	chet1 <- f22$y[cut(confidence[gte == 2], breaks=f22$x, labels=FALSE)]
##	chet2 <- f21$y[cut(confidence[gte == 2], breaks=f21$x, labels=FALSE)]

##	probWrongCall <- options@probWrongCall
##	term5 <- options@term5
##	gt.gte <- options@gt.gte
##	pTruthIsNormal[gte == 2] <- chet1*(1-term5["AB", "Normal"]) + chet2*term5["AB", "Normal"]
##	pTruthIsNormal[gte == 2] <- chet1*(1-gt.gte[2]) + chet2*gt.gte[2]

	##-------------------------------------------------------------------------
	##Calculate P(phat | Normal, HOM)
	##-------------------------------------------------------------------------	
	##In a normal region truth can be HOM or HET
	##P(phat | Normal, CHOM) is approx. P(phat | HET, CHOM) * P(HET|CHOM) + P(phat | HOM, CHOM) * P(HOM | CHOM)
##	chom1 <- f12$y[cut(confidence[gte == 1], breaks=f12$x, labels=FALSE)]
##	##chom3 <- 1-0.9999  ##P(HET|CHOM)
##	chom2 <- f11$y[cut(confidence[gte == 1], breaks=f11$x, labels=FALSE)]
##	##chom4 <- 0.9999    ##P(HOM|CHOM)
##	##probability that the true state is normal when genotype call is homozygous
##	pTruthIsNormal[gte == 1] <- chom1*(1-term5["AAorBB", "Normal"]) + chom2*term5["AAorBB", "Normal"]	
##	pTruthIsNormal[gte == 1] <- chom1*gt.gte[1] + chom2*(1-gt.gte[1])
	
	##-------------------------------------------------------------------------	
	##the stateProbability matrix is R x 2 (R = number of SNPs)
	##stateProbability is R x 2 (R rows and 2 columns)
	##gt.confidence.states <- options@gt.confidence.states
##	stateProbability <- stateProbability[, gt.confidence.states]
	
  	##-------------------------------------------------------------------------	
	##Calculate P(CHET| q) and P(CHOM| q), where q is
	##the true state.
  	##-------------------------------------------------------------------------		
	##For q = LOH, we make the following simplifying
	##assumptions:
	##1.  P(CHET | q=LOH) \approx P(CHET | HOM) 
	##2.  P(CHOM | q=LOH) \approx P(CHOM | HOM)
	##For q = Normal,
	##3.  P(CHET | q=Normal) \approx P(HET | q = Normal)
	##    estimate the latter from the data, or use a fixed number.  e.g., 0.3
	##4.  P(CHOM | q=Normal) \approx P(HOM | q = Normal)
	##    estimate the latter quantity from the data, or use a fixed number 0.7
	##P.CHET.LOH <- 1-P.CHOM.LOH
	##P.CHET.LOH <- 1-gte.state["L"]
	##P.CHET.Normal <- 1-P.CHOM.Normal
	##P.CHET.Normal <- 1-gte.state["N"]
	
	##*gte.state <- options@gte.state
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


getChromosomeArm <- function(snpset){
	data(chromosomeAnnotation, package="SNPchip", envir=environment())
	chrAnn <- as.matrix(chromosomeAnnotation)
	tmp <- position(snpset) <= chromosomeAnnotation[chromosome(snpset), "centromereStart"]
	chromosomeArm <- as.character(tmp)
	chromosomeArm[chromosomeArm == "TRUE"] <- "p"
	chromosomeArm[chromosomeArm == "FALSE"] <- "q"
	chromosomeArm
}



addFeatureData <- function(snpset){
	featureNames <- featureNames(snpset)

	if(!("arm" %in% fvarLabels(snpset))){
		arm <- getChromosomeArm(snpset)
	}
	fD <- data.frame(cbind(chromosome(snpset), arm=chromosomeArm), row.names=featureNames,
			 stringsAsFactors=FALSE)
	fD$position <- position(snpset) 
	colnames(fD) <- c("chromosome", "arm", "position")
	varMetadata <- data.frame(labelDescription=c("chromosome", "chromosomal arm", "physical position"), row.names=colnames(fD))
	featureData <- new("AnnotatedDataFrame", data=fD, varMetadata=varMetadata)
	featureData
}
	

##.calculateYlim <- function(object, op){
##  if("copyNumber" %in% ls(assayData(object))){
##    ##only print this if there is more than one sample or more than 1 chromosome to plot
##    if(length(unique(chromosome(object))) > 1 || ncol(object) > 1){
##      print("one.ylim is FALSE. Calculating ylim based on the percentiles of the copy number distribution")
##    }
##    if(op$log == "y"){
##      ylim <- range(copyNumber(object), na.rm=TRUE)
##    } else{
##      ##use 1 and 99 quantiles to avoid outliers
##      ylim <- c(quantile(copyNumber(object), prob=0.001, na.rm=TRUE),
##                quantile(copyNumber(object), prob=0.999, na.rm=TRUE))
##    }
##  } else{
##    ##ylimits for genotypes??
##    y <- .getY(object)
##    ##jitter the genotype calls
##    y <- jitter(y, amount=0.05)
##    ylim <- range(y)
##  }
##  ylim
##}




##For each of the breaks in xBreak, find the cytoband(s)
.locateCytoband <- function(breaks, cytoband){
  cytoband <- cytoband[cytoband[, "chrom"] == breaks["chromosome"], ]
  ##include cytobands that begin before the break and end in the break or after the break
  bands <- cytoband[cytoband[, "chromEnd"] > as.numeric(breaks["start"]) & cytoband[, "chromStart"] < as.numeric(breaks["last"]), "name"]
  bands <- paste(bands, collapse=", ")
}

##setMethod("calculateBreakpoints", "HmmParameter",
##          function(object, x, position, chromosome, states, digits=2, sampleNames, ...){
.calculatebreaks <- function(object, digits=2){
	if(class(object) != "HmmPredict") stop("arg object in .calcualtebreaks should be an object of class HmmPredict")
	position <- position(object)
	chromosome <- unique(chromosome(object))
	if(length(chromosome) > 1) stop("only 1 chromosome at a time")
	states <- states(object)
	sampleNames <- sampleNames(object)
	if(length(sampleNames) > 1) stop("only 1 sample at a time")
	x <- predictions(object)
	if(any(is.na(x))){
		Nmissing <- sum(is.na(x))
	} else Nmissing <- 0
	N <- length(x)
	colnames <- c("id", "chr", "state", "size", 
		      "N", "start", "last", "prev", "next")
	
	if(length(unique(x)) == 1){
		breaks <- data.frame(sampleNames,
				     chromosome,
				     states[unique(x)],
				     (max(position) - min(position))/1e6,
				     length(position),
				     min(position)/1e6,
				     max(position)/1e6,
				     NA,
				     NA)
		colnames(breaks) <- colnames
		return(breaks)
	}
	
	d <- diff(x)
	if(sum(abs(d), na.rm=TRUE) < 1){
		breaks <- data.frame(sampleNames,
				     chromosome,
				     states[unique(x)],
				     (max(position) - min(position))/1e6,
				     length(position),
				     min(position)/1e6,
				     max(position)/1e6,
				     NA, NA)
		colnames(breaks) <- colnames
		return(breaks)		    
	}
	index <- c(1, (2:N)[d != 0 & !is.na(d)])
	index <- cbind(index[-length(index)], index[2:length(index)])
	index[, 2] <- index[, 2] - 1
	lastBreak <- index[nrow(index), 2]
	index <- rbind(index, c(lastBreak+1, N))  

	##index is a matrix of indices in object where breaks occured.
	##Replace by physical position.
	physical.positions <- matrix(position[as.vector(index)], nrow(index), ncol(index))
	physical.positions <- physical.positions/1e6
	colnames(physical.positions) <- c("start", "last")
	size <- (physical.positions[, "last"] - physical.positions[, "start"])
	N <- index[, 2] - index[, 1] + 1
	physical.positions <- cbind(physical.positions, size, N)
	colnames(physical.positions)[3:4] <- c("size", "N")
	predictedState <- function(i, x){
		if(any(is.na(i))){
			stop("missing values")
		}
		pred <- x[(i[1]:i[2])]
		pred <- pred[!is.na(pred)]
		pred <- unique(pred)
		if(length(pred) > 1) {
			stop("predictions not unique")
		}
		pred
	}
	##The vector of predicted states (integer)
	ps <- t(apply(index, 1, predictedState, x))
	states <- states[as.numeric(ps)]
	breaks <- data.frame(physical.positions)
	breaks$state <- states
	breaks$chrom <- chromosome
	##breaks$id <- rep(sampleNames(object), nrow(breaks))
	
	##---------------------------------------------------------------------------
	##Find positions of adjacent SNPs
	##---------------------------------------------------------------------------	
	adjacentSnps <- function(x, position, chromosome, digits){
		start <- as.numeric(x["start"])*1e6
		last <- as.numeric(x["last"])*1e6
		chrom <- x["chrom"]
		if(any(position < start & chromosome == chrom)){
			prevSnp <- max(position[position < start & chromosome == chrom])
		} else prevSnp <- NA
		if(any(position > last & chromosome == chrom)){                 
			nextSnp <- min(position[position > last & chromosome == chrom])
		} else nextSnp <- NA
		if(!is.na(prevSnp)) prevSnp <- prevSnp/1e6
		if(!is.na(nextSnp)) nextSnp <- nextSnp/1e6
		adj <- c(prevSnp, nextSnp)
		adj
	}
	adjacent <- t(apply(breaks, 1, adjacentSnps, position=position,
			    chromosome=chromosome, digits=digits))
	colnames(adjacent) <- c("prev", "next")
	breaks <- cbind(breaks, adjacent)
	breaks$id <- rep(sampleNames, nrow(breaks))
	colnames(breaks) <- c("start", "last", "size", "N", "state", "chr", "prev", "next", "id")
	##Add columns for the name of the first SNP and the name of the last SNP
	breaks <- breaks[, colnames]
	breaks
}
	

calculateCnSE <- function(object,
                          referenceSet,
                          epsilon=0.1){
	if(min(copyNumber(object), na.rm=TRUE) > 0){
		print("Robust estimates of the standard error are on the log2 scale.")
		print("Transforming copy number in object to log2 scale")
		copyNumber(object) <- log2(copyNumber(object))
	} else{
		warning("Negative values in the copy number.  Assume that data in the copyNumber element of assayData has been suitably transformed and is approximately Gaussian")
	}
	if(!missing(referenceSet)){
		if(min(copyNumber(referenceSet), na.rm=TRUE) > 0){
			print("Transforming copy number in referenceSet to log2 scale")
			copyNumber(referenceSet) <- log2(copyNumber(referenceSet))
		} else{
			warning("Negative values in the copy number.  Calculations assume that data in the copyNumber element of assayData has been suitably transformed and is approximately Gaussian")
		}
	}
	if(missing(referenceSet)){
		print("using object to compute across sample standard deviations...assumes reasonable sample size")
		referenceSet <- object
	}
	require(genefilter) || stop("genefilter not available")
	is.autosome <- chromosome(object) %in% as.character(1:22) 
	object <- object[is.autosome, ]  
	referenceSet <- referenceSet[match(featureNames(object), featureNames(referenceSet)), ]
	robustSD <- function(X){
		diff(quantile(X, probs=c(0.16, (1-0.16)), na.rm=TRUE))/2 
	}
	within.sd <- apply(copyNumber(object), 2, robustSD)
	across.sd <- apply(copyNumber(referenceSet), 1, robustSD)
	##  across.sd <- rowSds(copyNumber(referenceSet), na.rm=TRUE)
	across.sd <- matrix(across.sd, nrow=nrow(object), ncol=ncol(object), byrow=FALSE)
	##scale across.sd by the median sd of the sample
	median.across.sd <- median(across.sd)
	std.across.sd <- across.sd/median.across.sd
	SE <- within.sd*std.across.sd
	SE[SE == 0] <- epsilon
	rownames(SE) <- featureNames(object)
	SE
}

scaleTransitionProbability <- function(object, SCALE=2, normalLabel="N"){##HmmOptions
	if(!(normalLabel %in% object@states)) stop("must supply the label for the normal (baseline) state")
	states <- object@states
	S <- length(states)
	if(S > 2){
		tau.scale <- matrix(1, S, S)
		rownames(tau.scale) <- colnames(tau.scale) <- states
		x <- seq(1, (S-1), by=0.01)
		y <- rep(NA, length(x))
		for(i in 1:length(x)){
			y[i] <- ((S-1)-x[i])/(S-2) 
		}
		z <- abs(x/y - SCALE)
		z <- z == min(z)
		x <- x[z]
		y <- y[z]
		tau.scale[upper.tri(tau.scale)] <- y
		tau.scale[lower.tri(tau.scale)] <- y
		i <- match(normalLabel, states)
		tau.scale[, i] <- x
		##by default, do not tau.scale the probabilities of leaving the normal state
		tau.scale[i, ] <- 1
		diag(tau.scale) <- 1
	} else{
		if(S == 2) tau.scale <- matrix(1, nrow=S, ncol=S)
	}
	return(tau.scale)
}

##Code this in C
viterbi <- function(beta, pi, tau, states, arm, tau.scale){
	S <- length(states)
	T <- length(tau) + 1
	##beta <- matrix(beta, nr=T, nc=S)
	delta <- psi <- matrix(NA, T, S)
	colnames(delta) <- colnames(psi) <- states
	delta[1, ] <- pi + beta[1, ]
	psi[1, ] <- rep(0, S)
	for(t in 2:T){
		if(arm[t] != arm[t-1]){
			##SNP t is on a different chromosome/chromosome arm
			delta[t, ] <- pi + beta[t, ]
			psi[t, ] <- rep(0, S)
			next()
		}
		AA <- matrix(tau[t-1], nr=S, nc=S)
		AA[upper.tri(AA)] <- AA[lower.tri(AA)] <- (1-tau[t-1])/(S-1)
		AA <- log(AA*tau.scale)
		for(j in 1:S){
			##Equation 105b
			delta[t, j] <- max(delta[t-1, ] + AA[, j]) + beta[t, j]
			psi[t, j] <- order(delta[t-1, ] + AA[, j], decreasing=TRUE)[1]
		}
	}
	Pstar <- max(delta[nrow(delta), ])
	qhat <- rep(NA, nrow(delta))
	T <- nrow(delta)
	names(qhat) <- rownames(delta)  
	qhat[T] <- order(delta[T, ], decreasing=TRUE)[1]
	for(t in (T-1):1){
		if(arm[t] != arm[t+1]){
			qhat[t] <- order(delta[t, ], decreasing=TRUE)[1]
		} else {
			qhat[t] <- psi[t+1, qhat[t+1]]
		}
	}
	qhat    
}

##this methods loops over snps for all samples, but is 15 times slower
.viterbi2 <- function(object, params){ ##HmmOptions, HmmParameters
	snpset <- object@snpset
	##	i <- match(featureNames(snpset(object)), dimnames(emission(params))[[1]])
	if(!(identical(featureNames(snpset), dimnames(emission(params))[[1]]))) stop("feature names should match")
	arm <- paste(chromosome(snpset), featureData(snpset)$arm, sep="")
	beta <- emission(params)
	S <- length(object@states)
	T <- nrow(snpset)
	
	tau <- genomicDistance(params)
	pi <- log(pi(params))
	tau.scale <- transitionScale(params)
	if(nrow(tau.scale) == 0){
		tau.scale <- matrix(1, S, S)
	}
	delta <- psi <- array(NA, dim=c(nrow(snpset), ncol(snpset), S))
	dimnames(delta) <- dimnames(psi) <- list(featureNames(snpset),
						 sampleNames(snpset),
						 states(object))
	pi <- matrix(pi, nrow=ncol(snpset), ncol=S, byrow=TRUE)
	delta[1, , ] <- pi + beta[1, , ]
	initial.psi <- matrix(0, nr=ncol(snpset), nc=S)
	psi[1, , ] <- initial.psi##rep(0, S)

	##handling of missing values for SNP t: the delta and psi
	##values for SNP t are assigned the delta and psi from the
	##previous SNP
	
	for(t in 2:T){
		if(t%%1000 == 0) cat(".")
		if(arm[t] != arm[t-1]){
			##SNP t is on a different chromosome/chromosome arm
			delta[t, , ] <- pi + beta[t, , ]
			psi[t, ,  ] <- initial.psi
			next()
		}
		AA <- matrix(tau[t-1], nr=S, nc=S)
		AA[upper.tri(AA)] <- AA[lower.tri(AA)] <- (1-tau[t-1])/(S-1)
		AA <- log(AA*tau.scale)
##		AA <- aperm(array(AA, dim=c(S, S, ncol(snpset))))
		##dimnames(AA) <- list(sampleNames(snpset), states(object), states(object))
		##Loop over states
		for(j in 1:S){
			##Equation 105b
			x <- matrix(delta[t-1, , ], nrow=ncol(snpset))
			y <- matrix(AA[, j], nrow=ncol(snpset), byrow=TRUE)
			z <- x + y
			delta[t, , j] <- apply(z, 1, "max") + beta[t, , j]
			f <- function(x) order(x, decreasing=TRUE)[1]
			psi[t, , j] <- apply(z, 1, f)
		}
##		if(any(is.na(delta[t, , ]))){
##			x <- matrix(delta[t, , ], nrow=dim(delta)[2])
##			samples <- which(apply(x, 1, function(x) any(is.na(x))))
##			delta[t, samples, ] <- delta[t-1, samples, ]
##			psi[t, samples, ] <- psi[t-1, samples, ]
##			next()
##		}
##		if(any(is.na(psi[t, , ]))){
##			x <- matrix(psi[t, , ], nrow=dim(psi)[2])
##			samples <- which(apply(x, 1, function(x) any(is.na(x))))
##			delta[t, samples, ] <- delta[t-1, samples, ]
##			psi[t, samples, ] <- psi[t-1, samples, ]
##		}
		if(t == T) cat("\n")
	}
	x <- matrix(delta[T, , ], nrow=dim(delta)[2])
	Pstar <- apply(x, 1, "max")
	qhat <- matrix(NA, nrow(delta), ncol(snpset))
##	rownames(qhat) <- dimnames(delta)[[1]]
##	colnames(qhat) <- sampleNames(snpset)
	qhat[T, ] <- apply(x, 1, function(x) order(x, decreasing=TRUE)[1])
	for(t in (T-1):1){
		if(arm[t] != arm[t+1]){
			x <- matrix(delta[t, , ], nrow=dim(delta)[2])
			qhat[t, ] <- apply(x, 1, function(x) order(x, decreasing=TRUE)[1])
		} else {
			qhat[t, ] <- psi[t+1, , qhat[t+1]]
		}
	}
	qhat
}

	

##How often is a SNP altered across samples
makeTable <- function(object, state, by, unit=1000, digits=3){
	if(missing(state)) stop("must specify state")
	if(missing(by)) by <- "chr"
	breaks <- breakpoints(object)
	X <- breaks[breaks[, "state"] == state, ]
	##split by chromosome
	X <- split(X, X[, by])
	##Frequency of alterations of a specific type per chromosome
	if(length(X) == 1){
		freq <- nrow(X[[1]])
	} else{
		freq <- sapply(X, nrow)
	}
	##Length of alterations
	med.L <- sapply(X, function(x) median(as.numeric(x[, "size"])*unit, na.rm=TRUE))
	avg.L <- sapply(X, function(x) mean(as.numeric(x[, "size"])*unit, na.rm=TRUE))
	sd.L <- sapply(X, function(x) sd(as.numeric(x[, "size"])*unit, na.rm=TRUE))
	##Number of SNPs per alteration
	med.N <- sapply(X, function(x) median(as.numeric(x[, "N"]), na.rm=TRUE))
	avg.N <- sapply(X, function(x) mean(as.numeric(x[, "N"]), na.rm=TRUE))
	sd.N <- sapply(X, function(x) sd(as.numeric(x[, "N"]), na.rm=TRUE))	
	stats <- as.data.frame(cbind(freq, med.L, avg.L, sd.L, med.N, avg.N, sd.N))
	colnames(stats) <- c("Freq", "med(length)", "avg(length)", 
			     "sd(length)", "med(n.snp)", "avg(n.snp)", 
			     "sd(n.snp)")
	stats[, c(3, 4, 6, 7)] <- round(stats[, c(3, 4, 6, 7)], digits)
	if("X" %in% rownames(stats)){
		rownames(stats)["X"] <- 23
	}
	stats <- stats[order(as.numeric(rownames(stats))), ]
	
	##Marginal stats
	X <- breaks[breaks[, "state"] == state, ]
	if(nrow(stats) <= 1){
		return(stats)
	}
	freq <- c(median(stats[, "Freq"]),
		  mean(stats[, "Freq"]),
		  sd(stats[, "Freq"]))
	med.L <- c(median(stats[, "med(length)"]),
		   mean(stats[, "med(length)"]),
		   sd(stats[, "med(length)"]))
	avg.L <- c(median(stats[, "avg(length)"]),
		   mean(stats[, "avg(length)"]),
		   sd(stats[, "avg(length)"]))
	avg.L <- round(avg.L, digits)
	sd.L <- c(median(stats[, "sd(length)"]),
		   mean(stats[, "sd(length)"]),
		   sd(stats[, "sd(length)"]))
	sd.L <- round(sd.L, digits)
	med.N <- c(median(stats[, "med(n.snp)"]),
		   mean(stats[, "med(n.snp)"]),
		   sd(stats[, "med(n.snp)"]))		
	avg.N <- c(median(stats[, "avg(n.snp)"]),
		   mean(stats[, "avg(n.snp)"]),
		   sd(stats[, "avg(n.snp)"]))
	avg.N <- round(avg.N, digits)
	sd.N <- c(median(stats[, "sd(n.snp)"]),
		   mean(stats[, "sd(n.snp)"]),
		   sd(stats[, "sd(n.snp)"]))
	sd.N <- round(sd.N, digits)
	overall <- cbind(freq, med.L, avg.L, sd.L, med.N, avg.N, sd.N)
	rownames(overall) <- c("median", "avg", "sd")
	stats <- rbind(stats, overall)
	stats
}








