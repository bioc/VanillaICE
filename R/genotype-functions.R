genotypeEmissionCrlmm <- function(object,
				  gt.conf,
				  cdfName,
				  prHetCalledHom,
				  prHetCalledHet,
				  prHomInNormal,
				  prHomInRoh){
	stopifnot(is(object, "matrix"))
	GT <- as.integer(object)
	##rm(object); gc()
	if(cdfName == "pd.genomewidesnp.6"){
		annotation <- "genomewidesnp6"
	} else annotation <- cdfName
	loader(paste(annotation, "Conf.rda", sep=""), .vanillaIcePkgEnv, "VanillaICE")
	hapmapP <- getVarInEnv("reference")
##	prHetCalledHom <- hmm.params[["prHetCalledHom"]]
##	prHetCalledHet <- hmm.params[["prHetCalledHet"]]
##	prHomInNormal <- hmm.params[["prHomInNormal"]]
##	prHomInRoh <- hmm.params[["prHomInRoh"]]
	if(length(cdfName) < 1) stop("must specify annotation")
	##data(list=paste(annotation, "Conf", sep=""), package="VanillaICE", envir=environment())
	if(length(prHomInNormal) == nrow(gt.conf)){  ##convert to vector
		prHomInNormal <- as.numeric(matrix(prHomInNormal, nrow(gt.conf), ncol(gt.conf), byrow=FALSE))
	} else prHomInNormal <-  rep(prHomInNormal, length(GT))
	hapmapP[, 2] <- 1-exp(-hapmapP[, 2]/1000)
	##p = 1-exp(-X/1000)
	##1000*log(1-p)=X
	##confidence <- 1-exp(-GTconf/1000)
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
	##GT <- as.integer(genotypes)
	##confidence <- as.numeric(confidence)
	##confidence <- as.numeric(gt.conf)
	confidence <- gt.conf
	rm(gt.conf); gc()
	pTruthIsNormal <- pTruthIsRoh <- rep(NA, length(GT))
	confidence[confidence==0] <- 0.01 ##Otherwise, NA's result
	hom <- which(GT == 1 | GT == 3)
	observedPcalledHom <- cut(confidence[hom], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[hom] <- f11$y[observedPcalledHom]
	het <- which(GT == 2)
	observedPcalledHet <- cut(confidence[het], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[het] <- f21$y[observedPcalledHet]
	##-------------------------------------------------------------------------
	##Calculate P(phat | Normal, HOM)
	##-------------------------------------------------------------------------
	chet1 <- f22$y[cut(confidence[het], breaks=f22$x, labels=FALSE)]
	chet2 <- f21$y[cut(confidence[het], breaks=f21$x, labels=FALSE)]
	##term5[1]=P(true genotype is HET | genotype call is AB, state is normal)
	pTruthIsNormal[het] <- chet1*prHetCalledHet + chet2*(1-prHetCalledHet)
	##chom1=called homozygous truth heterozygous
	chom1 <- f12$y[cut(confidence[hom], breaks=f12$x, labels=FALSE)]
	##chom2=called homozygous truth homozygous
	chom2 <- f11$y[cut(confidence[hom], breaks=f11$x, labels=FALSE)]
	##chom4 <- 0.9999    ##P(HOM|CHOM)
	##probability that the true state is HOM when genotype call is homozygous
	##prHetCalledHom = P(true genotype is HET | calls is AA or BB, state is normal)
	pTruthIsNormal[hom] <- chom1*prHetCalledHom + chom2*(1-prHetCalledHom)
	fNormal <- fLoh <- rep(NA, length(GT))
	fNormal[hom] <- prHomInNormal[hom] * pTruthIsNormal[hom]
	fNormal[het] <- (1-prHomInNormal[het]) * pTruthIsNormal[het]
	fLoh[hom] <- prHomInRoh * pTruthIsRoh[hom]
	fLoh[het] <- (1-prHomInRoh) * pTruthIsRoh[het]
	f <- array(NA, dim=c(nrow(object), ncol(object), 2)) ##dimnames=list(featureNames(object),
							##     sampleNames(object),
							  ##   c("normal", "ROH")))
	f[, , 1] <- matrix(fNormal, length(object), ncol(object))
	f[, , 2] <- matrix(fLoh, length(object), ncol(object))
	dimnames(f)[[3]] <- c("normal", "ROH")
	f[f  == 0] <- min(f[f > 0], na.rm=TRUE)
	f <- log(f)
	return(f)
}

computeLoglikFromViterbi2 <- function(object, chrom, id, pos, log.it=TRUE){
	viterbiSequence <- viterbiStatePath(object)
	S <- nStates(object)
	rl <- Rle(viterbiSequence)
	starts <- start(rl)
	LLR <- rep(NA,  length(starts))
	if(log.it){
		log.emission <- matrix(log(emission(object)), object@numberFeatures, nStates(object))
		log.initial <- log(initialProb(object))
	} else {
		log.emission <- matrix(emission(object), object@numberFeatures, nStates(object))
		log.initial <- initialProb(object)
	}
	p <- c(NA, transitionProb(object))
	##p <- c(NA, as.numeric(viterbiResults[["tau"]]))
	c1 <- normal2altered(object)
	c2 <- altered2normal(object)
	c3 <- altered2altered(object)
	lP.N2N <- log(1-((1-p)*(S-1)*c1)) ##probability normal -> normal
	lP.N2A <- log((1-p)*c1) ##probability normal -> altered
	P.A2A <- sapply(1-((1-p)*(c2+(S-2)*c3)), function(x) max(x, 0.01))
	lP.A2A <- log(P.A2A) ## probability altered to same altered state
	lP.A2N <- log((1-p)*c2) ##probability altered -> normal
	lP.A2Astar <- log((1-p)*c3) ## probability altered -> different altered state
	for(k in seq_along(starts)){
		LLR[k] <- computeLoglikForRange(from=start(rl)[k],
						to=end(rl)[k],
						viterbiSequence=viterbiSequence,
						normalIndex=normalStateIndex(object),
						log.initial=log.initial,
						log.emission=log.emission,
						lP.N2N=lP.N2N,
						lP.N2A=lP.N2A,
						lP.A2N=lP.A2N,
						lP.A2A=lP.A2A)
	}
	start.index <- start(rl)
	end.index <- end(rl)
	start <- pos[start.index]
	end <- pos[end.index]
	numMarkers <- width(rl)
	states <- viterbiSequence[start.index]
	ir <- IRanges(start=start, end=end)
	charChrom <- oligoClasses::integer2chromosome(chrom)
	##sl <- getSequenceLengths(genomeBuild(object))[charChrom]
	rd <- GRanges(seqnames=paste("chr", charChrom, sep=""),
		      ir,
		      state=as.integer(states),
		      numberProbes=numMarkers,
		      LLR=LLR)
##	names(rd) <- id
##	rd <- RangedDataHMM(ir,
##			    chromosome=rep(chrom, length(LLR)),
##			    sampleId=rep(id, length(LLR)),
##			    state=states,
##			    coverage=numMarkers,
##			    startIndexInChromosome=start.index,
##			    endIndexInChromosome=end.index,
##			    LLR=LLR)
	return(rd)
}

gtEmissionFromMatrix <- function(object,
				 S,
				 gt.conf,
				 is.snp,
				 normalIndex=3L,
				 rohIndex=c(2L, 4L),
				 ICE=FALSE,
				 prGtHom=c(0.7, 0.99),
				 prGtMis=rep(1/S,S),
				 prHetCalledHom=0.001,
				 prHetCalledHet=0.995,
				 prHomInNormal=0.8,
				 prHomInRoh=0.999,
				 annotationPkg,
				 log.it=TRUE,
				 ...){
	if(!ICE){
		p <- prGtHom
		prGenotypeMissing <- prGtMis
		stopifnot(length(p) == S)
		if(!is.numeric(object)) stop("genotypes must be integers (1=AA, 2=AB, 3=BB) or NA (missing)")
		emission <- array(NA, dim=c(nrow(object), ncol(object), S))
		missingGT <- any(is.na(object))
		for(s in seq_len(S)){
			tmp <- object
			tmp[tmp == 1 | tmp == 3] <- p[s]
			tmp[tmp == 2] <- 1-p[s]
			index1 <- is.na(tmp) & !is.snp
			index2 <- is.na(tmp) & is.snp
			if(missingGT){
				tmp[index2] <- prGenotypeMissing[s]
				tmp[index1] <- 1/S
			}
			emission[, , s] <- tmp
		}
		if(log.it) emit <- log(emission) else emit <- emission
	} else {
		not.valid <- invalidGtConfidence(gt.conf)
		if(any(not.valid)){
			stop("Invalid genotype confidence scores.\n",
			     "\tIf ICE is TRUE, all confidence scores must be between 0 and 1")
		}
		logemit <- array(NA, dim=c(nrow(object), ncol(object), S))
		pkg <- strsplit(annotationPkg, "Crlmm")[[1]][[1]]
		tmp <- genotypeEmissionCrlmm(object, gt.conf=gt.conf,
					     cdfName=pkg,
					     prHetCalledHom,
					     prHetCalledHet,
					     prHomInNormal,
					     prHomInRoh)
		logemit[, , rohIndex] <- tmp[, , "ROH"]
		logemit[, , -rohIndex] <- tmp[, , "normal"]
		if(log.it) emit <- logemit else emit <- exp(logemit)
	}
	return(emit)
}

invalidGtConfidence <- function(x){
	is.na(x) | x < 0 | x > 1 | is.nan(x) | is.infinite(x)
}

viterbiForSnpSet <- function(gt, S=2L,
			     is.snp, pos, gt.conf,
			     chrom,
			     ICE=FALSE,
			     rohIndex=2L,
			     annotationPkg,
			     prGtHom,
			     prGtMis=rep(1/S,S),
			     prHetCalledHom=0.001,
			     prHetCalledHet=0.995,
			     prHomInNormal=0.8,
			     prHomInRoh=0.999,
			     normalIndex=1L,
			     TAUP=1e8,
			     ...){ ## additional arguments passed to Viterbi (e.g., TAUP, normal2altered...)
	log.beta.gt <- gtEmissionFromMatrix(object=gt,
					    S=S,
					    is.snp=is.snp,
					    pos=pos,
					    gt.conf=gt.conf,
					    ICE=ICE,
					    annotationPkg=annotationPkg,
					    prGtHom=prGtHom,
					    prGtMis=prGtMis,
					    prHetCalledHom=prHetCalledHom,
					    prHetCalledHet=prHetCalledHet,
					    prHomInNormal=prHomInNormal,
					    prHomInRoh=prHomInRoh,
					    rohIndex=2L)
	transitionProb <- computeTransitionProb(pos, S=S, TAUP=TAUP)
	## construct ViterbiSet
	J <- ncol(gt)
	viterbiList <- vector("list", J)
	for(j in seq_len(J)){
		tmp <- new("Viterbi",
			   numberFeatures=nrow(gt),
			   numberStates=S,
			   emission=as.matrix(as.numeric(log.beta.gt[, j, ])),
			   normalIndex=normalIndex,
			   transitionProb=transitionProb,
			   ...)
		viterbiList[[j]] <- fit(tmp)
	}
	id <- object <- NULL
	rdlist <- foreach(object = viterbiList,
			  id=colnames(gt)) %do% {
				  computeLoglikFromViterbi2(object=object,
							    chrom=chrom,
							    id=id,
							    pos=pos,
							    log.it=FALSE)
			  }
	res <- GRangesList(rdlist)
	names(res) <- colnames(gt)
	return(res)
}

setMethod("fit", signature(object="Viterbi"), function(object){
	fitViterbi(object)
})

fitViterbi <- function(object){
	tmp <- .C("viterbi",
		  emission(object),
		  initialProb(object),
		  transitionProb(object), ## log scale?
		  object@arm,
		  nStates(object),
		  object@numberFeatures,
		  viterbiStatePath(object),
		  delta(object),
		  object@normal2altered,
		  object@altered2normal,
		  object@altered2altered,
		  normalStateIndex(object),
		  pAA(object))
	new("Viterbi",
	    emission=emission(object),
	    initialProb=initialProb(object),
	    transitionProb=transitionProb(object),
	    arm=object@arm,
	    numberStates=nStates(object),
	    numberFeatures=object@numberFeatures,
	    viterbiStatePath=tmp[[7]],
	    delta=tmp[[8]],
	    normal2altered=object@normal2altered,
	    altered2normal=object@altered2normal,
	    altered2altered=object@altered2altered,
	    normalIndex=normalStateIndex(object),
	    pAA=pAA(object))
}

computeLoglikForRange <- function(from, to,
				  viterbiSequence, normalIndex,
				  log.initial,
				  log.emission,
				  lP.A2N,
				  lP.N2A,
				  lP.N2N,
				  lP.A2A){
	index <- seq(from, to)
	thisState <- unique(viterbiSequence[index])
	if(thisState == normalIndex) return(0)
	first.index <- min(index)
	last.index <- max(index)
	T <- length(viterbiSequence)
	## 6 Rules  (1 < t < T)
	## 1.  index=1
	## 2.  index=t
	## 3.  index=t,t+1
	## 4.  index=T
	## 5.  index=1,2
	## 6,  index=T-1, T
	##------
	## 1. index=1
	if(first.index == last.index & last.index==1){
		##note the last term cancels
		logLik.vit <- log.initial[thisState]+log.emission[1, thisState] + lP.A2N[2] + log.emission[2, normalIndex]
		logLik.null <- log.initial[normalIndex]+log.emission[1, normalIndex] + lP.N2N[2] + log.emission[2, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	##2 index = t
	if(length(index) == 1 & first.index > 1 & last.index < T){
		##note the last term cancels
		logLik.vit <- sum(lP.N2A[index] + log.emission[index, thisState]) + lP.A2N[last.index+1]+log.emission[last.index+1, normalIndex]
		logLik.null <- sum(lP.N2N[index] + log.emission[index, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	##if T=t+1?
	## 3: index = t, ..., t*  , t>1, t* < T, t* != t
	if(first.index != last.index & first.index > 1 & last.index < T){
		index2 <- index[-1]
		logLik.vit <- lP.N2A[first.index] +
			sum(lP.A2A[index2]) +
				lP.A2N[last.index+1] +
					log.emission[first.index, thisState] +
						sum(log.emission[index2, thisState]) +
							log.emission[last.index+1, normalIndex]
		logLik.null <-
			sum(lP.N2N[index]) +
				lP.N2N[last.index+1]  +
					sum(log.emission[index, normalIndex]) +
						log.emission[last.index+1, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	## 4: index = T
	if(first.index == last.index & last.index == T){
		logLik.vit <- lP.N2A[T] + log.emission[T, thisState]
		logLik.null <- lP.N2N[T] + log.emission[T, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	## 5: index = 1, 2, ...
	if(first.index != last.index & first.index == 1 & last.index < T){
		index2 <- index[-1]## t=2, ...., t*
		logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
		logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	if(first.index != last.index & first.index == 1 & last.index == T){
		index2 <- index[-1]## t=2, ...., t*
		logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
		logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	## 6: index = t, ...T
	if(first.index != last.index & last.index == T){
		index2 <- index[-1]
		logLik.vit <- lP.N2A[first.index] + log.emission[first.index, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
		logLik.null <- lP.N2N[first.index] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
		LLR <- logLik.vit - logLik.null
		return(LLR)
	}
	stop("none of conditions in likRatio function satisfied")
}

checkAnnotationForICE <- function(object){
	if(!annotation(object) %in% icePlatforms()){
		stop("ICE is TRUE, but hapmap crlmm confidence scores for ", annotation(object), " are not available.")
	}
}


icePlatforms <- function(){
	c("pd.genomewidesnp.6",
	  "genomewidesnp6",
	  "genomewidesnp6Crlmm",
	  "pd.mapping250k.nsp",
	  "pd.mapping250k.sty",
	  "pd.mapping250k.nsp, pd.mapping250k.sty")
}

setMethod("gtEmission", signature(object="matrix"),
	  function(object, gt.conf,
		   is.snp,
		   rohIndex=c(2L, 4L),
		   S=6L,
		   annotationPkg,
		   ICE=FALSE,
		   prGtHom=c(0.7,0.99),
		   prGtMis=rep(1/S, S),
		   prHetCalledHom=0.001,
		   prHetCalledHet=0.995,
		   prHomInNormal=0.8,
		   prHomInRoh=0.999, ...){
		  gtEmissionFromMatrix(object, gt.conf=gt.conf,
				       is.snp=is.snp,
				       rohIndex=rohIndex,
				       S=S,
				       annotationPkg=annotationPkg, ICE=ICE,
				       prGtHom=prGtHom,
				       prGtMis=prGtMis,
				       prHetCalledHom=prHetCalledHom,
				       prHetCalledHet=prHetCalledHet,
				       prHomInNormal=prHomInNormal,
				       prHomInRoh=prHomInRoh,...)
	  })

gtEmissionFromMatrix <- function(object,
				 S,
				 gt.conf,
				 is.snp,
				 normalIndex=3L,
				 rohIndex=c(2L, 4L),
				 ICE=FALSE,
				 prGtHom=c(0.7, 0.99),
				 prGtMis=rep(1/S,S),
				 prHetCalledHom=0.001,
				 prHetCalledHet=0.995,
				 prHomInNormal=0.8,
				 prHomInRoh=0.999,
				 annotationPkg,
				 log.it=TRUE,
				 ...){
	if(!ICE){
		p <- prGtHom
		prGenotypeMissing <- prGtMis
		stopifnot(length(p) == S)
		if(!is.numeric(object)) stop("genotypes must be integers (1=AA, 2=AB, 3=BB) or NA (missing)")
		emission <- array(NA, dim=c(nrow(object), ncol(object), S))
		missingGT <- any(is.na(object))
		for(s in seq_len(S)){
			tmp <- object
			tmp[tmp == 1 | tmp == 3] <- p[s]
			tmp[tmp == 2] <- 1-p[s]
			index1 <- is.na(tmp) & !is.snp
			index2 <- is.na(tmp) & is.snp
			if(missingGT){
				tmp[index2] <- prGenotypeMissing[s]
				tmp[index1] <- 1/S
			}
			emission[, , s] <- tmp
		}
		if(log.it) emit <- log(emission) else emit <- emission
	} else {
		not.valid <- invalidGtConfidence(gt.conf)
		if(any(not.valid)){
			stop("Invalid genotype confidence scores.\n",
			     "\tIf ICE is TRUE, all confidence scores must be between 0 and 1")
		}
		logemit <- array(NA, dim=c(nrow(object), ncol(object), S))
		pkg <- strsplit(annotationPkg, "Crlmm")[[1]][[1]]
		tmp <- genotypeEmissionCrlmm(object, gt.conf=gt.conf,
					     cdfName=pkg,
					     prHetCalledHom,
					     prHetCalledHet,
					     prHomInNormal,
					     prHomInRoh)
		logemit[, , rohIndex] <- tmp[, , "ROH"]
		logemit[, , -rohIndex] <- tmp[, , "normal"]
		if(log.it) emit <- logemit else emit <- exp(logemit)
	}
	return(emit)
}
