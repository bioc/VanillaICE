rowMAD <- function(x, y, ...){
	notna <- !is.na(x)
	sds <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(sds)
}

robustSds <- function(x, takeLog=FALSE, ...){
	if(!is.matrix(x)) stop("x is not a matrix")
	if(takeLog) x <- log2(x)
##	if(ncol(x) > 3){
##		sds1 <- rowMAD(x, ...)
##		sds1 <- matrix(sds1, nrow(x), ncol(x))
##		sds2 <- apply(x, 2, "mad", ...)
##		df <- ncol(x)
##		##sds2 <- sds2/median(sds2, na.rm=TRUE)
##		##sds2 <- sds2/min(sds2, na.rm=TRUE)
##		##sds <- t(t(sds1) * sds2)
##	} else {
	sds <- apply(x, 2, "mad", ...)
	sds <- matrix(sds, nrow(x), ncol(x), byrow=TRUE)
	dimnames(sds) <- dimnames(x)
	return(sds)
}

viterbi.wrapper <- function(log.emission,
			    log.initial,
			    transitionPr,
			    arm,
			    S,
			    T,
			    result,
			    delta,
			    normal2altered=1,
			    altered2normal=1,
			    altered2altered=1,
			    normalIndex,
			    pAA){
	tmp <- list(as.matrix(as.double(as.matrix(exp(log.emission)))),##beta
		    as.double(as.matrix(exp(log.initial))),##initialP
		    as.matrix(as.double(transitionPr)),
		    as.integer(arm),##arm
		    as.integer(S),##number of states
		    as.integer(T),##number of loci
		    result,##placeholder for results
		    as.matrix(as.double(delta)),##delta
		    normal2altered=normal2altered,##c1
		    altered2normal=altered2normal,##c2
		    altered2altered=altered2altered,##c3
		    as.integer(normalIndex),
		    as.double(pAA))##normalState
	res <- .C("viterbi",
		  log.emission=tmp[[1]],
		  log.initial=tmp[[2]],
		  tau=tmp[[3]],
		  arm=tmp[[4]],
		  S=tmp[[5]],
		  T=tmp[[6]],
		  viterbiSeq=tmp[[7]],
		  delta=tmp[[8]],
		  N2A=tmp[[9]],
		  A2N=tmp[[10]],
		  A2A=tmp[[11]],
		  normalIndex=tmp[[12]],
		  pAA=tmp[[13]])
	return(res)
}

.getArm <- function(chrom, pos, genome){
	if(is.integer(chrom)) chrom <- paste("chr", integer2chromosome(chrom), sep="")
	path.gap <- system.file("extdata", package="SNPchip")
	gaps <- readRDS(list.files(path.gap, pattern=paste("gap_", genome, ".rda", sep=""), full.names=TRUE))
	centromere.starts <- start(gaps)
	names(centromere.starts) <- seqnames(gaps)
	centromere.starts <- centromere.starts[chrom]
	arm <- ifelse(pos <= centromere.starts, "p", "q")
	paste(chrom, arm, sep="")
}

getChromosomeArm <- function(object){
	chrom <- chromosome(object)
	pos <- position(object)
	res <- .getArm(chrom, pos, genomeBuild(object))
	return(res)
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

computeLoglik <- function(viterbiResults,
			  log.initial,
			  log.emission=log.emission,
			  states,
			  normalIndex,
			  nNotMissing,
			  c1, c2, c3,
			  physical.position,
			  CHR,
			  sample.name){
	viterbiSequence <- viterbiResults[["viterbiSeq"]]
	S <- length(states)
	rl <- Rle(viterbiSequence)
	starts <- start(rl)
	LLR <- rep(NA,  length(starts))
	log.emission <- matrix(viterbiResults[["log.emission"]], nNotMissing, S)
	##** The NA is to stagger the transition probabilities by 1
	##  -- this way, the same index can be used to multiply the transition and emission probabilities
	p <- c(NA, as.numeric(viterbiResults[["tau"]]))
	lP.N2N <- log(1-((1-p)*(S-1)*c1)) ##probability normal -> normal
	lP.N2A <- log((1-p)*c1) ##probability normal -> altered
	P.A2A <- sapply(1-((1-p)*(c2+(S-2)*c3)), function(x) max(x, 0.01))
	lP.A2A <- log(P.A2A) ## probability altered to same altered state
	lP.A2N <- log((1-p)*c2) ##probability altered -> normal
	lP.A2Astar <- log((1-p)*c3) ## probability altered -> different altered state
	##For each seqment, compute the likelihood ratio
	for(k in seq_along(starts)){
		LLR[k] <- computeLoglikForRange(from=start(rl)[k],
						to=end(rl)[k],
						viterbiSequence=viterbiSequence,
						normalIndex=normalIndex,
						log.initial=log.initial,
						log.emission=log.emission,
						lP.N2N=lP.N2N,
						lP.N2A=lP.N2A,
						lP.A2N=lP.A2N,
						lP.A2A=lP.A2A)
	}
	start.index <- start(rl)
	end.index <- end(rl)
	##pos <- position(object)[I]
	##this is tricky since we've added an index to force a segment for each arm.
	start <- physical.position[start.index]
	end <- physical.position[end.index]
	##numMarkers <- unlist(numMarkers)
	numMarkers <- width(rl)
	states <- viterbiSequence[start.index]
	ir <- IRanges(start=start, end=end)
	rangedData <- RangedData(ranges=ir,
				 chromosome=rep(CHR, length(LLR)),
				 sampleId=rep(sample.name, length(LLR)),
				 state=states,
				 coverage=numMarkers,
				 LLR=LLR)
	return(rangedData)
}



##iterbi <- function(##object,
##		    log.initial,
##		    normal2altered,
##		    altered2normal,
##		    altered2altered,
##		    normalIndex,
##		    states,
##		    log.E,
##		    TAUP, ...){
##	sns <- colnames(log.E)
##	if(is.null(sns)) stop("no dimnames for log.emission")
##	##log.initial <- hmm.params$log.initialPr
##	##verbose <- hmm.params$verbose
###	normal2altered <- hmm.params$n2a
###	altered2normal <- hmm.params$a2n
###	altered2altered <- hmm.params$a2a
##	if(normal2altered <= 0) stop("normal2altered must be > 0")
##	if(altered2normal <= 0) stop("altered2normal must be > 0")
##	if(altered2altered <= 0) stop("altered2altered must be > 0")
###	##
###	## arm can contain NA's for invalid chromosomes or NA's
###	arm <- getChromosomeArm(object)
###	arm <- as.integer(as.factor(arm))
##	##normalIndex <- hmm.params$normalIndex##[["normalIndex"]]
##	if(normalIndex < 1 | normalIndex > dim(log.E)[[3]]){
##		stop("normalIndex in hmm.params not valid")
##	}
###	TAUP <- hmm.params$tau
##	c1 <- normal2altered
##	c2 <- altered2normal
##	c3 <- altered2altered
##	##TT <- nrow(object)
##	index <- which(!is.na(arm))
##	arm <- arm[index]
##	TT <- length(index)
##	log.E <- log.E[index, , , drop=FALSE]
##	##all.equal(log.E[1:10, 1, ], LE[1:10, 1, ])
##	object <- object[index, ]
##	##
##	##states <- hmm.params$states
##	names(log.initial) <- as.character(states)
##	S <- length(states)
##	delta <- matrix(as.double(0), nrow=TT, ncol=S)
##	rangedData <- vector("list", ncol(log.E))
##	for(j in 1:ncol(log.E)){
##		rD <- vector("list", length(unique(arm)))
##		missingE <- rowSums(is.na(log.E[, j, ])) > 0
##		notFinite <- rowSums(!is.finite(log.E[, j, ])) > 0
##		missingE <- missingE | notFinite
##		for(a in seq_along(unique(arm))){
##			logicalNotMissing <- I <- arm == a  & !missingE
##			if(sum(logicalNotMissing) < 2) next()
##			nNotMissing <- T <- sum(logicalNotMissing)
##			physical.position <- position(object)[logicalNotMissing]
##			transitionPr <- exp(-2 * diff(physical.position)/TAUP)
##			##is the lower bound a function of normal2altered, altered2normal, altered2altered?
##			minimum <- 1-1/((S-1)*c1) + 0.01
##			transitionPr <- pmax(transitionPr, minimum)
##			if(any(transitionPr < 0 | transitionPr > 1)) stop("Transition probabilities not in [0,1].  Order object by chromosome and physical position")
##			result <- rep(as.integer(0), T)
##			viterbiResults <- viterbi.wrapper(log.emission=log.E[I, j, ],
##							  log.initial=log.initial,
##							  transitionPr=transitionPr,
##							  arm=arm[I],
##							  S=S,
##							  T=T,
##							  result=result,
##							  delta=delta[I, ],
##							  normal2altered=normal2altered,
##							  altered2normal=altered2normal,
##							  altered2altered=altered2altered,
##							  normalIndex=normalIndex,
##							  pAA=rep(0, S^2))
##			rD[[a]] <- computeLoglik(viterbiResults,
##						 log.initial=log.initial,
##						 log.emission=log.E[I, j, ],
##						 states=states,
##						 normalIndex=normalIndex,
##						 nNotMissing=nNotMissing,
##						 c1=c1, c2=c2, c3=c3,
##						 physical.position=physical.position,
##						 CHR=unique(chromosome(object)[logicalNotMissing]),
##						 sample.name=sns[j])
##		}
##		notnull <- !sapply(rD, is.null)
##		rD <- rD[notnull]
###		if(length(rD)==1) {
###			rangedData[[j]] <- rD[[1]]
###		} else{
###			rangedData[[j]] <- stack(RangedDataList(rD))
###		}
##		rD <- GRangesList(rD)
##	}
##	if(length(rangedData) == 1) {
##		rangedData <- rangedData[[1]]
##	} else {
##		rangedData <- stack(RangedDataList(rangedData))
##	}
##	return(rangedData)
##

stackRangedData <- function(object){
	if(is(object, "list")){
		if(length(object)==1){
			object <- object[[1]]
		} else {
			object <- RangedDataList(object)
			object <- stack(object)
			ix <- match("sample", colnames(object))
			if(length(ix) > 0) object <- object[, -ix]
		}
		if(is(object, "RangedDataHMM")) return(object)
	}
	rangedData <- RangedDataHMM(ranges=ranges(object),
				    chromosome=object$chrom,
				    sampleId=object$sampleId,
				    state=object$state,
				    coverage=object$coverage,
				    LLR=object$LLR)
	return(rangedData)
}
##
##stackRangedDataHMM <- function(object){
##
##
##}

rbaf <- function(genotypes, sigma, epsilon, states){
	baf <- matrix(NA, nrow(genotypes), ncol(genotypes))
	Ns <- table(genotypes)
	a <- pnorm(0, mean=0, sd=sigma)
	b <- pnorm(1, mean=0, sd=sigma)
	I <- runif(Ns[1], 0, 1) > epsilon
	baf[genotypes==1] <- I*qnorm(a+runif(Ns[1], 0, b-a), mean=0, sd=sigma) + (1-I)*runif(Ns[1], 0, 1)
	I <- runif(Ns[2], 0, 1) > epsilon
	baf[genotypes==2] <- I*rnorm(Ns[2], mean=0.5, sd=sigma*2) + (1-I)*runif(Ns[2], 0, 1)
	a <- pnorm(0, mean=1, sd=sigma)
	b <- pnorm(1, mean=1, sd=sigma)
	I <- runif(Ns[3], 0, 1) > epsilon
	baf[genotypes==3] <- I*qnorm(a+runif(Ns[3], 0, b-a), mean=1, sd=sigma) + (1-I) * runif(Ns[3], 0, 1)

	## assume 1/2 are hets to make it easy
	ndup <- sum(states==5)
	ndup.het <- ceiling(ndup/2)
	ndup.hom <- ndup-ndup.het

	index5 <- which(states==5)
	index25 <- sample(index5, ceiling(ndup.het/2))
	index5 <- setdiff(index5, index25)

	index75 <- setdiff(index5, floor(ndup.het/2))
	indexhom <- setdiff(index5, index75)

	n25 <- length(index25)
	n75 <- length(index75)
	nhom <- length(index5)
	I <- runif(n25, 0, 1) > epsilon
	baf[index25] <- I*rnorm(n25, mean=1/3, sd=sigma*2) + (1-I)*runif(n25, 0, 1)
	I <- runif(n75, 0, 1) > epsilon
	baf[index75] <- I*rnorm(n75, mean=2/3, sd=sigma*2) + (1-I)*runif(n75, 0, 1)
	a <- pnorm(0, mean=1, sd=sigma)
	b <- pnorm(1, mean=1, sd=sigma)
	baf[indexhom] <- qnorm(a+runif(Ns[3], 0, b-a), mean=1, sd=sigma)
	rownames(baf) <- rownames(genotypes)
	return(baf)
}


centerAutosomesAt <- function(x, at, ...){
	stopifnot(!missing(at))
	marker.index <- which(chromosome(x) <= 23)
	cn <- copyNumber(x)[marker.index, , drop=FALSE]
	meds <- apply(cn, 2, median, na.rm=TRUE)
	cn.cen <- sweep(cn, 2, meds) + at
	copyNumber(x)[marker.index, ] <- cn.cen
	return(x)
}



##hmm.setup <- function(object, ...){  ## whether the save the emission probabilities
##	res <- HmmOptionList(object=object, ...)
##	res <- as(res, "HmmOptionList")
##	return(res)
##}
hmm.setup <- function(...) .Deprecated("hmm.setup function is deprecated.  hmm will run directly on a class of oligoSnpSet, SnpSet, CopyNumberSet, or BeadStudioSet. See ?hmm for details.")


setMethod("update", "environment", function(object, ...){
	if(length(ls(object)) == 0) stop("nothing to update")
	hmm(object=object[["locusset"]],
	    states=object[["states"]],
	    mu=object[["mu"]],
	    probs=object[["probs"]],
	    takeLog=object[["takeLog"]],
	    initialP=object[["initialP"]],
	    returnSegments=object[["returnSegments"]],
	    TAUP=object[["TAUP"]],
	    verbose=object[["verbose"]],
	    ice=object[["ice"]],
	    envir=object)
})


invalidCnConfidence <- function(x){
	is.na(x) | x <= 0 | is.nan(x) | is.infinite(x)
}

invalidGtConfidence <- function(x){
	is.na(x) | x < 0 | x > 1 | is.nan(x) | is.infinite(x)
}

getSds <- function(object, na.rm=TRUE){
	cn.conf <- cnConfidence(object)
	chrom <- chromosome(object)
	stopifnot(all(chrom <= 24))
	notvalid <- invalidCnConfidence(cn.conf)
	CN <- copyNumber(object)
	if(any(notvalid)){
		sds <- .getSds(CN, chrom)
		##if(verbose) message("cnConfidence missing.  Using MAD")
##		marker.index <- which(chrom < 23)
##		if(length(marker.index) == 0){
##			sds <- matrix(NA, nrow(cn.conf), ncol(cn.conf))
##			## sex chromosomes
##			marker.index.list <- split(seq_len(nrow(CN)), chrom)
##			for(i in seq_along(marker.index.list)){
##				marker.index <- marker.index.list[[1]]
##				tmp <- apply(CN[marker.index, , drop=FALSE], 2, mad, na.rm=TRUE)
##				sds[marker.index, ] <- matrix(tmp, length(marker.index), ncol(CN), byrow=TRUE)
##			}
##		}  else {  ## autosomes present
##			s <- apply(CN[marker.index, , drop=FALSE], 2, mad, na.rm=TRUE)
##			sds <- matrix(s, nrow(CN), ncol(CN), byrow=TRUE)
##		}
		valid <- !notvalid
		if(any(valid)){
			sds[valid] <- 1/cn.conf[valid]
		}
	} else { ## all valid confidence scores
		sds <- 1/cn.conf
	}
	notvalid <- invalidCnConfidence(sds)
	stopifnot(any(!notvalid))
	return(sds)
}

.getSds <- function(CN, chrom){
	nr <- nrow(CN)
	nc <- ncol(CN)
	if(!missing(chrom)){
		marker.index <- which(chrom < 23)
		if(length(marker.index) > 0){
			CN <- CN[marker.index, , drop=FALSE]
			s <- apply(CN, 2, mad, na.rm=TRUE)
			sds <- matrix(s, nrow(CN), ncol(CN), byrow=TRUE)
		} else {
			sds <- matrix(NA, nr, nc)
			## sex chromosomes
			marker.index.list <- split(seq_len(nrow(CN)), chrom)
			for(i in seq_along(marker.index.list)){
				marker.index <- marker.index.list[[1]]
				tmp <- apply(CN[marker.index, , drop=FALSE], 2, mad, na.rm=TRUE)
				sds[marker.index, ] <- matrix(tmp, length(marker.index), ncol(CN), byrow=TRUE)
			}
		}
	} else { ## chromosome missing
		sds <- apply(CN, 2, mad, na.rm=TRUE)
		##sds <- matrix(s, nr, nc, byrow=TRUE)
	}
	return(sds)
}


validChromosomeIndex <- function(object){
	index <- which(chromosome(object) <= 24 & !is.na(chromosome(object)) & !is.na(position(object)))
	if(length(index) < 1){
		stop("Chromosome must be 1-24 and physical position \n
                      can not be missing.  See chromosome() and position().\n
                      See integer2chromosome(1:24) for integer codes used by\n
                      VanillaICE.")
	}
	return(index)
}

probabilityOutlier <- function(cn, k=3, sigmas, verbose=TRUE){
	## outlier ~ N(0, sigma1), cn ~ N(0, sigma2), sigma2 << sigma1
	## lik= prod_i=1^N Pr(outlier) N(0, sigma1) + (1-Pr(outlier)) N(0, sigma2)
	if(length(cn) > k){
		rmeds <- runmed(cn, k)
		delta <- cn-rmeds
	} else delta <- as.numeric(cn-median(cn, na.rm=TRUE))
	if(missing(sigmas)) {
		sd.delta <- mad(delta, na.rm=TRUE)
		if(sd.delta < 0.01){
			sd.delta <- sd(delta, na.rm=TRUE)
			if(sd.delta < 0.01){
				stop("standard deviation of delta is near 0 in probabilityOutlier function. Probably too few markers.")
			}
		}
		sigmas <- c(sd.delta*3, sd.delta)  ## set variance for outlier component to be 3 times the normal component
	}
	mu <- 0
	tau <- 0.01
	epsilon <- 2; counter <- 1
	while(epsilon > 0.01){
		## two component Gaussian mixture
		## sigma[1] is the variance for the outlier component
		gamma <- (tau * dnorm(delta, mu, sigmas[1]))/ (tau * dnorm(delta, mu, sigmas[1]) + (1-tau)*dnorm(delta, mu, sigmas[2]))
		## gamma near 1 is likely an outlier
		## gamma near 0 is likely not an outlier
		tau.next <- mean(gamma)
		epsilon <- abs(tau.next - tau)
		tau <- tau.next
		counter <- counter+1
		if(counter > 10) break()
	}
	## plot(delta, gamma)
	return(gamma)
}

##probabilityOutlierBaf <- function(b, centers, sigmas, lower, upper, tau=rep(1/length(centers), length(centers))){
##	## outlier ~ N(0, sigma1), cn ~ N(0, sigma2), sigma2 << sigma1
##	## lik= prod_i=1^N Pr(outlier) N(0, sigma1) + (1-Pr(outlier)) N(0, sigma2)
##	epsilon <- 2; counter <- 1
##	while(epsilon > 0.01){
##		## two component Gaussian mixture
##		## sigma[1] is the variance for the outlier component
##		##gamma <- (tau * dnorm(delta, mu, sigmas[1]))/ (tau * dnorm(delta, mu, sigmas[1]) + (1-tau)*dnorm(delta, mu, sigmas[2]))
##		num <- tau[1]*tnorm(b, centers[1], sigmas[1],
##				    lower=lower[1], upper=upper[1])
##		tmp <- vector("list", length(centers))
##		for(i in seq_along(centers)){
##			tmp[[i]] <- tau[i]*tnorm(b, centers[i],
##						 sigmas[i],
##						 lower=lower[i],
##						 upper=upper[i])
##		}
##		gamma <- vector("list", length(centers))
##		for(i in seq_along(centers[-1])){
##			gamma[[i]] <- tmp[[1]]/sum(tmp[-1])
##		## gamma near 1 is likely an outlier
##		## gamma near 0 is likely not an outlier
##		tau.next <- mean(gamma)
##		epsilon <- abs(tau.next - tau)
##		tau <- tau.next
##		counter <- counter+1
##		if(counter > 10) break()
##	}
##	## plot(delta, gamma)
##	return(gamma)
##}




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




####
#### update mu and sigma using the forward backward variables
####
#### Should updateMu drop copy-neutral ROH (as is currently done?)
####
##updateMu <- function(x,
##		     mu,
##		     sigma,
##		     normalIndex,
##		     rohIndex=normalIndex+1,
##		     is.log,
##		     limits,
##		     fv,
##		     bv,
##		     constrainMu){
##	if(nUpdates==0) return(mu)
##	min.sd <- mad(x, na.rm=TRUE)/2
##	sigma <- pmax(sigma, min.sd)
##	s <- sigma
##	L <- S <- length(mu)
####	if(length(dup.index) > 0){
####		mu <- unique(mu)
####		L <- length(mu)
####	} else L <- S
##	## fix normal copy number
##	mu[normalIndex] <- median(x, na.rm=TRUE)
##	mu[rohIndex] <- mu[normalIndex]
##	p <- matrix(0, L, 2)
##	p[, 1] <- 0.99
##	p[, 2] <- 1-p[, 1]
##	##sigma <- rep(sigma, L)
##	g2 <- g1 <- matrix(NA, length(x), L)
##	missing.fv <- missing(fv)
##	if(!missing.fv){
##		h <- matrix(fv*bv, length(x), length(mu))
##		m <- rowSums(h, na.rm=TRUE)
##		h <- h/m
##		rm(m, fv, bv);
##	}
##	d1 <- matrix(NA, length(x), L)
##	na.index <- which(is.na(x))
##	anyNA <- length(na.index) > 0
##	nr <- length(x)
##
##	d.norm <- dnorm(x, mean=matrix(mu, nr, L, byrow=TRUE), sd=matrix(sigma, nr, L, byrow=TRUE))
##	d.norm <- t(t(d.norm) * p[, 1])
##	d.unif <- dunif(x, matrix(limits[1], nr, L), matrix(limits[2], nr, L))
##	d.unif <- t(t(d.unif) * p[, 2])
##	d1 <- d.norm/(d.norm+d.unif)
##	if(length(na.index) > 0)
##		d1[na.index, ] <- dunif(rep(0, length(na.index)), limits[1], limits[2])
##	if(missing.fv){
##		g1 <- d1
##		g2 <- 1-d1
##	} else {
##		g1 <- h * d1
##		g2 <- h * (1-d1)
##	}
##	##
##	## While the outlier component of the mixture is the
##	## same for each state, the probability of being an
##	## outlier will differ.  This makes sense as a small
##	## value should have a relatively large outlier
##	## probability for the normal distribution and
##	## relatively small for a deletion
##	##
##	total.g1 <- apply(g1, 2, sum, na.rm=TRUE)
##	## denom.eq52:
##	denom.eq52 <- apply(g1 + g2, 2, sum, na.rm=TRUE)
##	## p[i, j] is the expected number of times in state i
##	## using mixture component k relative to the expected
##	## number of times in state i
##	##
##	p.new <- p
##	p.new[, 1] <- total.g1/denom.eq52
##	p.new[, 2] <- 1-p.new[,1]
##	sigma.new <- mu.new <- rep(NA, L)
##	I <- seq_len(L)
##	tmp <- g1 * x
##	mu.new <- colSums(tmp, na.rm=TRUE)/total.g1
##	mu.new[-rohIndex] <- constrainMu(mu.new[-rohIndex], is.log)
##	mu.new[rohIndex] <- mu.new[normalIndex]
##	mu.new <- makeNonDecreasing(mu.new)
##	xx <- matrix(x, nr, L)
##	tmp <- g1*t(t(xx) - mu.new)^2
##	sigma.new <- sqrt(colSums(tmp, na.rm=TRUE)/total.g1)
##	sigma.new <- pmax(sigma.new, min.sd)
##	sigma.new <- constrainSd(sigma.new)
##	mu <- mu.new
##	sigma <- sigma.new
##	p <- p.new
##	return(list(mu, sigma, p))
##}


tnorm <- function(x, mean, sd, lower=0, upper=1){
       phi <- function(x, mu, sigma) dnorm(x, mu, sigma)
       ## cdf of standard normal
       Phi <- function(x, mu, sigma) pnorm(x, mu, sigma)
       res <- phi(x, mean, sd)/(Phi(upper, mean, sd)-Phi(lower, mean, sd))
##       ind <- which(x < lower | x > upper)
##       if(any(ind)){
##               res[ind] <- 0
##       }
       res
}

trnorm <- function(x){
	function(mu, sigma, upper=1, lower=0){
		dnorm(x, mu, sigma)/(pnorm(upper, mu, sigma) - pnorm(lower, mu, sigma))
	}
}

##trnorm <- function(x, j){
##	isMatrix <- is.matrix(x)
##	if(isMatrix) x <- x[, j]
##	##if(!isMatrix) x <- as(x, "matrix")
##	##nr <- nrow(x); nc <- ncol(x)
##	function(mu, sigma, upper=1, lower=0){
##		##Mu <- matrix(mu, nr, nc, byrow=TRUE)
##		##Sigma <- matrix(sigma, nr, nc, byrow=TRUE)
##		dnorm(x, Mu, Sigma)/(pnorm(upper, Mu, Sigma) - pnorm(lower, Mu, Sigma))
##	}
##}

trnorm2 <- function(x, j){
	isMatrix <- is.matrix(x)
	if(!isMatrix)
		x <- as(x, "matrix")
	nr <- nrow(x); nc <- ncol(x)
	x <- t(x)
	function(mu, sigma, upper=1, lower=0){
		Mu <- matrix(mu, nc, nr, byrow=FALSE)
		Sigma <- matrix(sigma, nc, nr, byrow=FALSE)
		dnorm(x, Mu, Sigma)/(pnorm(upper, Mu, Sigma) - pnorm(lower, Mu, Sigma))
	}
}


getTau <- function(TAUP, pos, S){
	taus <- exp(-2*diff(pos)/TAUP)
	tau.min <- 1-1/((S-1)) + 0.01
	taus <- pmax(taus, tau.min)
}





viterbi3 <- function(arm, pos, chrom, LE, log.initial,
		     states, TAUP=1e8, normalIndex=3, id){
	stopifnot(is(LE, "matrix"))
	S <- length(states)
	na.LE <- is.na(LE)
	if(any(na.LE)){
		missingE <- rowSums(na.LE) > 0
		notFinite <- rowSums(!is.finite(LE[, ])) > 0
		missingE <- missingE | notFinite
		I <- !missingE
		pos <- pos[I]
		arm <- arm[I]
		LE <- LE[I, ]
		if(is(arm, "character"))
			arm <- as.integer(as.factor(arm))
	}
	T <- nrow(LE)
	qhat <- rep(0L, T)
	delta <- matrix(as.double(0), nrow=T, ncol=S)
	taus <- getTau(TAUP, pos, S)
	viterbiResults <- viterbi.wrapper(log.emission=LE,
					  log.initial=log.initial,
					  transitionPr=taus,
					  arm=as.integer(as.factor(arm)),
					  S=S,
					  T=T,
					  result=qhat,
					  delta=delta,
					  normal2altered=1,
					  altered2normal=1,
					  altered2altered=1,
					  normalIndex=normalIndex,
					  pAA=rep(0, S^2))
	rd <- computeLoglik(viterbiResults,
			    log.initial=log.initial,
			    log.emission=LE,
			    states=states,
			    normalIndex=normalIndex,
			    nNotMissing=nrow(LE),
			    c1=1, c2=1, c3=1,
			    physical.position=pos,
			    CHR=unique(chrom),
			    sample.name=id)
	return(rd)
}

getColClasses <- function(filename, lrr.colname, baf.colname){
	tmp <- read.table(filename,
			  row.names=NULL,
			  header=TRUE,
			  stringsAsFactors=FALSE,
			  sep="\t", nrows=50)
	j <- grep(lrr.colname, colnames(tmp))
	k <- grep(baf.colname, colnames(tmp))
	if(length(j) == 0 || length(k)==0)
		stop("lrr.colname or baf.colname not in header")
	colClasses <- as.character(sapply(tmp[1, ], class))
	index <- setdiff(seq_along(colClasses), c(1, j, k))
	colClasses[index] <- rep("NULL", length(index)) ## don't read in the other columns
	colClasses
}

## copied from MinimumDistance
read.bsfiles <- function(path="", filenames, ext="", row.names=1,
			 sep="\t",
			 lrr.colname="Log.R.Ratio",
			 baf.colname="B.Allele",
			 drop=FALSE,
			 colClasses,
			 nrows=1.8e6,
			 ...){
	if(path != ""){
		fnames <- file.path(path, paste(filenames, ext, sep=""))
	} else fnames <- paste(filenames, ext, sep="")
	stopifnot(all(file.exists(fnames)))
	if(missing(colClasses)){
		colClasses <- getColClasses(fnames[1], lrr.colname, baf.colname)
	}
	for(i in seq_along(fnames)){
		##cat(".")
		tmp <- read.table(fnames[i],
				  row.names=row.names,
				  sep=sep,
				  nrows=nrows,
				  header=TRUE,
				  stringsAsFactors=FALSE,
				  colClasses=colClasses,
				  check.names=FALSE,
				  comment.char="", ...)
		tmp <- as.matrix(tmp)
		if(i==1){
			dat <- array(NA, dim=c(nrow(tmp), 2, length(fnames)))
			j <- grep(lrr.colname, colnames(tmp))
			k <- grep(baf.colname, colnames(tmp))
			stopifnot(length(j)==1)
			stopifnot(length(k)==1)
			if(!drop){
				dimnames(dat) <- list(rownames(tmp),
						      c("lrr", "baf"),
						      basename(filenames))
			}
			##lrr.data <- matrix(NA, nrow(tmp), length(filenames))
			##baf.data <- matrix(NA, nrow(tmp), length(filenames))
		}
		dat[, 1, i] <- tmp[, j]
		dat[, 2, i] <- tmp[, k]
	}
	##cat("\n")
	return(dat)
}

artificialData <- function(states, nmarkers){
	data(oligoSetExample, package="oligoClasses")
	oligoSet <- chromosomePositionOrder(oligoSet)
	pos <- position(oligoSet)
	state.path <- rep(states, nmarkers)
	copynumber <- rep(2, length(state.path))
	copynumber[state.path==2] <- 1.5##bias
	copynumber[state.path==5] <- 2.5##bias
	genotypes <- rep(NA, length(copynumber))
	gt <- rmultinom(n=length(copynumber), size=1, prob=rep(1/3,3))
	genotypes[gt[1, ] == 1] <- 1L
	genotypes[gt[2, ] == 1] <- 2L
	genotypes[gt[3, ] == 1] <- 3L
	genotypes[state.path==4 | state.path==2] <- 1L
	genotypes <- as.matrix(genotypes)
	## make signal fairly obvious
	sigmas <- rgamma(length(copynumber), 4, scale=0.05)
	b <- rbaf(as.matrix(genotypes), sigma=0.01, epsilon=0.001, states=state.path)
	dat <- as.matrix(rnorm(length(state.path), mean=copynumber, sd=sigmas))
	i <- seq_along(state.path)
	rownames(dat) <- rownames(genotypes) <- featureNames(oligoSet)[i]
	##pos <- seq(1, by=3e3, length.out=length(copynumber))
	object <- new("oligoSnpSet",
		      copyNumber=integerMatrix(dat, scale=100),
		      call=as.matrix(genotypes),
		      callProbability=snpCallProbability(oligoSet)[i, , drop=FALSE])
	##baf(object) <- b
	assayDataElement(object, "baf") <- integerMatrix(b, scale=1000)
	##df <- data.frame(position=pos, chromosome=rep(1L, length(pos)), isSnp=
	fData(object)$position <- as.integer(pos[i])
	fData(object)$chromosome <- 1L
	fData(object)$isSnp <- TRUE
	return(object)
}

mask <- function(query, subject){
	##index <- matchMatrix(findOverlaps(query, subject))[, 1]
	index <- queryHits(findOverlaps(query, subject))
	query <- query[-index, ]
}


getProbB <- function(cdfname, featurenames){
	cdfpath <- system.file("extdata", package=cdfname)
	if(file.exists(file.path(cdfpath, "pb_gw6.rda"))){
		load(file.path(cdfpath, "pb_gw6.rda"))
		probB <- rep(NA, length(featurenames))
		pb <- pb[names(pb) %in% featurenames]
		index <- match(names(pb), featurenames)
		stopifnot(!any(is.na(index)))
		probB[index] <- pb
	} else probB <- rep(0.5, length(featurenames))
	return(probB)
}

keyOffFirstFile <- function(filename, cdfname, genome, lrr.colname, baf.colname, ...){
	## read in one file
	## return feature matrix in chromosome, position order
	dat <- read.bsfiles(filenames=filename, lrr.colname=lrr.colname, baf.colname=baf.colname, ...)
	cdfpath <- system.file("extdata", package=cdfname)
	load(file.path(cdfpath, paste("snpProbes_", genome, ".rda", sep="")))
	load(file.path(cdfpath, paste("cnProbes_", genome, ".rda", sep="")))
	snpProbes <- get("snpProbes")
	cnProbes <- get("cnProbes")
	features <- rbind(snpProbes, cnProbes)
	keep.index <- which(rownames(features) %in% rownames(dat))
	features <- features[keep.index, ]

	index.order <- order(features[, "chrom"], features[, "position"])
	features <- features[index.order, ]

	issnp <- as.logical(rownames(features) %in% rownames(snpProbes))
	probB <- as.integer(getProbB(cdfname, rownames(features))*100)
	arm <- .getArm(features[, "chrom"], features[, "position"], genome)
	index <- match(rownames(features), rownames(dat))

	identical(rownames(dat)[index], rownames(features))
	##features2 <- cbind(features, issnp, probB, arm, index)
	features2 <- data.frame(chrom=features[, "chrom"],
				position=features[, "position"],
				probB=probB,
				isSnp=issnp,
				arm=arm,
				index=index)
	##colnames(features2) <- c(colnames(features), "isSnp", "probB", "arm", "index")
	return(features2)
}

## ad-hoc.  Do we want to put priors on the means?
constrainMu <- function(mu, is.log){
	if(is.log){
		is.ratio <- all.equal(mu[3], 0, tolerance=0.2) == TRUE
		if(is.ratio){
			mu[3] <- ifelse(mu[3] < -0.2, -0.2, mu[3])
			mu[3] <- ifelse(mu[3] > 0.2, 0.2, mu[3])
			mu[1] <- ifelse(mu[1] > -1, -1, mu[1])
			mu[2] <- ifelse(mu[2] > -0.25, -0.25, mu[2])
			mu[4] <- ifelse(mu[4] < 0.25, 0.25, mu[4])
			mu[4] <- ifelse(mu[4] > 0.6, 0.6, mu[4])
			mu[5] <- ifelse(mu[5] < 0.65, 0.65, mu[5])
		} else {
			mu[3] <- ifelse(mu[3] < 0.75, 0.75, mu[3])
			mu[3] <- ifelse(mu[3] > 1.25, 1.25, mu[3])
			mu[1] <- ifelse(mu[1] > -0.5, -0.5, mu[1])
			mu[2] <- ifelse(mu[2] > 0.7, 0.7, mu[2])
			mu[4] <- ifelse(mu[4] < 1.25, 1.25, mu[4])
			mu[4] <- ifelse(mu[4] > 2, 2, mu[4])
			mu[5] <- ifelse(mu[5] < 2, 2, mu[5])
		}
	} else {
		mu[1] <- ifelse(mu[1] > 0.8, 0.8, mu[1])
		mu[2] <- ifelse(mu[2] > 1.7, 1.7, mu[2])
		mu[2] <- ifelse(mu[2] < 0.9, 0.9, mu[2])
		mu[3] <- ifelse(mu[3] > 2.2, 2.2, mu[3])
		mu[3] <- ifelse(mu[3] < 1.8, 1.8, mu[3])
		mu[4] <- ifelse(mu[4] < 2.3, 2.3, mu[4])
		mu[4] <- ifelse(mu[4] > 3, 3, mu[4])
		mu[5] <- ifelse(mu[5] < 3, 3, mu[5])
	}
	return(mu)
}

constrainSd <- function(sigma){
	sigma[1] <- min(3*sigma[3], max(sigma[1], sigma[3]))
	sigma[6] <- min(3*sigma[3], max(sigma[6], sigma[3]))
	sigma[5] <- min(3*sigma[3], max(sigma[5], sigma[3]))
	sigma[2:4] <- sigma[3]
	return(sigma)
}


computeTransitionProb <- function(x, TAUP, S, tauMAX=1-5e-6){
	p <- exp(-2*diff(x)/TAUP)
	minimum <- 1-1/((S-1)) + 0.01
	p <- pmax(p, minimum)
	p <- pmin(p, tauMAX)
	return(as.matrix(p))
}


copyNumberLimits <- function(is.log){
	if(is.log) return(c(-2.5, 3)) else return(c(0,10))
}

thresholdCopyNumber <- function(object, limits){
	object <- pmin(object, limits[2])
	object <- pmax(object, limits[1])
	object
}


centerCopyNumber <- function(object, is.snp){
	snp.index <- which(is.snp)
	mu.snp <- apply(object[snp.index, , drop=FALSE], 2, median, na.rm=TRUE)
	object[snp.index, ] <- sweep(object[snp.index, , drop=FALSE],  2, mu.snp)
	if(any(!is.snp)){
		np.index <- which(!is.snp)
		mu.np <- apply(object[np.index, , drop=FALSE], 2, median, na.rm=TRUE)
		object[np.index, ] <- sweep(object[np.index, , drop=FALSE], 2, mu.np)
	}
	object
}

guessCnScale <- function(x){
	is.log <- if(all(x >= 0, na.rm=TRUE)) FALSE else TRUE
	is.log
}

logCnStates <- function(normalCn, is.ratio) {
	if(missing(is.ratio))
		is.ratio <- all.equal(normalCn, 0, tolerance=0.3) == TRUE
	if(is.ratio){
		return(c(-2, -0.5, 0, 0, 0.4, 0.8))
	} else{
		is.absolute <- all.equal(normalCn, 1, tolerance=0.2)==TRUE
		if(is.absolute){
			return(c(-1.5, 0, 1, 1, 1.58, 2))
		} else{
			stop("if copy number is on the log base 2 scale, the copy number for the normal state should be either 0 or 1 depending on whether the copy number is relative or absolute, respectively.")
		}
	  }
}

copyNumberStates <- function(normalCn) {
	is.ratio <- all.equal(normalCn, 1, tolerance=0.2) == TRUE
	if(is.ratio){
		return(c(0, 1/2, 1, 1, 3/2, 4/2))
	} else{
		is.absolute <- all.equal(normalCn, 2, tolerance=0.2)==TRUE
		if(is.absolute){
			return(c(0, 1, 2, 2, 3, 4))
		} else{
			stop("median copy number is not near 1 or 2")
		}
	  }
}

rescale <- function(x, l, u){
	b <- 1/(u-l)
	a <- l*b
	(x+a)/b
}

makeNonDecreasing <- function(x){
	d <- diff(x)
	if(all(d >= 0)) return(x)
	index <- which(d < 0)
	if(index[1] == 1){
		x[1] <- x[2]
		index <- index[-1]
		if(length(index) ==0)
			return(x)
	}
	l <- length(index)
	if(l == length(x)-1){
		i <- index[l]
		x[length(x)] <- x[length(x)-1]
		index <- index[-l]
		if(length(index) == 0) return(x)
	}
	x[index+1] <- x[index]
	return(x)
}

stackGRangesList <- function(fit, build){
	L <- length(fit[[1]])
	gr <- list()
	for(i in seq_len(L)){
		rdl <- lapply(fit, "[[", i)
		tmp <- stack(GRangesList(rdl))
		gr[[i]] <- tmp[, -grep("sample", colnames(values(tmp)))]
	}
	sl <- getSequenceLengths(build)
	rd <- GRangesList(gr)#, seqlengths=sl)
	nms <- unique(as.character(runValue(unlist(seqnames(rd)))))
	seqlengths(rd) <- sl[match(nms, names(sl))]
	names(rd) <- names(fit[[1]])
	metadata(rd) <- list(genome=build)
	rd
}


initialSigma <- function(sds, J, S){
	sigmas <- matrix(sds, J, S, byrow=FALSE)
	## homozygous deletions often have a higher variance.
	sigmas[, 1] <- sigmas[, 1] * 4
	sigmas[, 2] <- sigmas[, 2] * 2
	sigmas
}

##initialParams <- function(J){
##	nms <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
##	mus <- function(J) {
##		x <- matrix(c(0, 1/4, 1/3, 0.5, 2/3, 3/4, 1), J, 7, byrow=TRUE)
##		colnames(x) <- nms
##		x
##	}
##	sds <- function(J){
##		x <- matrix(c(0.02, rep(0.05, 5), 0.02), J, 7, byrow=TRUE)
##		colnames(x) <- nms
##		x
##	}
##	p <- function(J){
##		p <- replicate(J, matrix(c(0.99, 0.01), S, 2, byrow=TRUE), simplify=FALSE)
##	}
##	list(mu=mus, sd=sds, p=p)
##}

generatorFun <- function(r, b, gt, is.snp, cnStates,
			 normalIndex, TAUP, limits, center,
			 prOutlierBAF, p.hom, position, is.log, computeLLR, chrom){
	S <- length(cnStates)
	nc <- ncol(r)
	nr <- nrow(r)
	np.index <- which(!is.snp)
	names(prOutlierBAF)[1] <- "prOutlier"
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")

	## params for copy number
	sds <- apply(r, 2, mad, na.rm=TRUE)
	if(center) r <- centerCopyNumber(r, is.snp) + cnStates[normalIndex]
	r <- thresholdCopyNumber(r, limits=limits)

	initialCnParams <- function(j){
		mus <- cnStates
		sigmas <- rep(sds[j], S)
		##p <- matrix(c(0.99,  0.01), S, 2, byrow=TRUE)
		p <- 0.01
		paramsCN <- list(mu=mus, sigmas=sigmas, p=p)
	}

	## params for BAF
	allele.prob <- getPrB()
	initialBafParams <- function(){
		musBAF <- c(0, 1/4, 1/3, 0.5, 2/3, 3/4, 1)
		sdsBAF <- c(0.02, rep(0.05, 5), 0.02)
		names(musBAF) <- names(sdsBAF) <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
		paramsBAF <- list(mus=musBAF, sigmas=sdsBAF, prOutlier=prOutlierBAF)
	}

	g2 <- g1 <- matrix(NA, nr, S)
	d1 <- matrix(NA, nr, S)
	isHom <- b < 0.02 | b > 0.98
	indexHom <- list()
	for(j in seq_len(nc)) indexHom[[j]] <- which(isHom[, j])
	emitr <- emitb <- matrix(NA, nr, S)
	##r.isna <- is.na(r)
	##r.hasna <- any(r.isna)
	cn.dunif <- dunif(r, limits[1], limits[2])
	Gtotal <- matrix(NA, nr, 4)
	updateBafEmission <- function(bf, params, j){
		mus <- params[["mus"]]
		sds <- params[["sigmas"]]
		##prOutlier <- params$prOutlier$initial
		prOutlier <- params$prOutlier$prOutlier
		computeEmitPr <- function(dlist, prOutlier, allele.prob){
			sumd <- 0
			q <- 1-prOutlier
			for(i in seq_along(dlist)){
				sumd <- sumd + dlist[[i]] * q*allele.prob[i]
			}
			##p <- matrix(prOutlier, nc, nr, byrow=TRUE)
			t(sumd+prOutlier)
		}
		trNormal <- trnorm(bf)

		##---------------------------------------------------------------------------
		## copy number 1
		dA <- trNormal(mus["A"], sds["A"])
		dB <- trNormal(mus["B"], sds["B"])
		emitb[, 2] <- computeEmitPr(list(dA,dB), prOutlier, allele.prob[, c("A", "B")])
		emitb[, 4] <- emitb[, 2]

		##---------------------------------------------------------------------------
		## copy number 2
		dAB <- trNormal(mus["AB"], sds["AB"])
		emitb[, 3] <- computeEmitPr(list(dA, dAB, dB), prOutlier=prOutlier, allele.prob=allele.prob[, c("AA", "AB", "BB")])

		##---------------------------------------------------------------------------
		## copy number 3
		dAAB <- trNormal(mus["AAB"], sds["AAB"])
		dABB <- trNormal(mus["ABB"], sds["ABB"])
		emitb[, 5] <- computeEmitPr(list(dA, dAAB, dABB, dB), prOutlier, allele.prob[, c("AAA", "AAB", "ABB", "BBB")])

		##---------------------------------------------------------------------------
		## copy number 4
		dAAAB <- trNormal(mus["AAAB"], sds["AAAB"])
		dABBB <- trNormal(mus["ABBB"], sds["ABBB"])
		emitb[, 6] <- computeEmitPr(list(dA, dAAAB, dAB, dABBB, dB), prOutlier, allele.prob[, c("AAAA", "AAAB", "AABB", "ABBB", "BBBB")])

		##
		## small p.hom makes homozygous genotypes less informative
		##
		if(p.hom < 1){
			i <- indexHom[[j]]
			if(length(i) > 0) emitb[i, c(3, 5, 6)] <- (1-p.hom)*emitb[i, 2] + p.hom*emitb[i, c(3, 5, 6)]
		}
		if(length(np.index) > 0) emitb[np.index, ] <- 1
		nas <- anyMissing(emitb)
		if(nas) emitb[is.na(emitb)] <- 1
		##index <- which(rowSums(is.infinite(emitb)) > 0)
		##if(any(is.infinite(emitb))) browser()
		return(emitb)
	}

	## might try dt here instead
	mydnorm <- function(x){
		function(mean, sd){
			mus <- matrix(mean, nr, S, byrow=TRUE)
			sds <- matrix(sd, nr, S, byrow=TRUE)
			dnorm(x, mus, sds)
		}
	}
	updateCnEmission <- function(cn, params, j){
		mus <- params[["mu"]]
		sds <- params[["sigmas"]]
		p <- params[["p"]]
		q <- 1-p
		normalDens <- mydnorm(cn)
		d <- normalDens(mus, sds)
		emitr <- q*d + p*cn.dunif[, j]
		emitr.na <- is.na(emitr)
		hasna <- any(emitr.na)
		if(hasna) emitr[emitr.na] <- 1
		return(emitr)
	}

	theoreticalMeans <- c(0,1/4, 1/3, 1/2, 2/3, 3/4, 1)
	names(theoreticalMeans) <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
	updateBafParams <- function(bf, params, fv, bv, j) {
		##hasna <- anyMissing(bf)
		##if(hasna){
		weightedMean <- function(w, sumWeights) sum(w*bf, na.rm=TRUE)/sumWeights
		weightedSd <- function(w, sumWeights, mu) sqrt((sum(w*(bf-mu)^2, na.rm=TRUE))/sumWeights)
		##		} else {
		##			weightedMean <- function(w, sumWeights) sum(w*bf)/sumWeights
		##			weightedSd <- function(w, sumWeights, mu) sqrt((sum(w*(bf-mu)^2))/sumWeights)
		##		}
		mus <- params[["mus"]]
		sigmas <- params[["sigmas"]]
		pOut <- params[["prOutlier"]]
		pOutlierMax <- pOut[["max"]]
		pOutlierMaxROH <- pOut[["maxROH"]]
		p <- pOut[["prOutlier"]]
		p <- pmin(p, pOutlierMax)
		##p[2] <- min(p[2], pOutlierMaxROH)
		p <- 1-p ## probability not an outlier

		## forward / backward variables
		h <- matrix(fv*bv, nr, S)
		m <- rowSums(h, na.rm=TRUE)
		h <- h/m
		calculateGamma <- function(d, h, allele.prob){
			sumd <- 0
			## conditional on not an outlier
			G <- length(d)
			## integrate over G possible genotypes
			## (marginal probability for CN given not an outlier)
			for(j in seq_len(G)){
				d[[j]] <- allele.prob[j]*d[[j]]
				sumd <- sumd+d[[j]]
			}
			## integrate out the outlier indicator to get the marginal probability for the copy number
			## (marginal probability for CN)
			##denom <- p*sumd + (1-p)
			lapply(d, function(d) d*h/sumd)
		}
		updateBMoments <- function(gammas, nms){
			## denom = q.out * [pr(AA) * phi(baf, muAA, sigmaAA) + pAB*phi(baf, muAB, sigmaAB) + pBB*phi(baf, muBB, sigmaBB)) + p.out
			totalWeights <- sapply(gammas, sum, na.rm=TRUE)
			isZero <- totalWeights == 0
			mus <- mapply(weightedMean, w=gammas, sumWeights=totalWeights)
			## for many near-zero bafs, it seems that the weight for genotype 'AA' is near zero...
			sigmas <- mapply(weightedSd, w=gammas, sumWeights=totalWeights, mu=mus)
			if(any(isZero)){
				## then the baf corresponding to the genotype is probably not observed
				## the mean should be sampled from the prior (or simply set to the theoretical value)
				gt <- nms[isZero]
				mus[isZero] <- theoreticalMeans[gt]
				sigmas[isZero] <- params$sigmas[gt]
			}
			moments <- cbind(mus, sigmas)
			rownames(moments) <- nms
			moments
		}
		nms <- names(mus)
		trNormal <- trnorm(bf)

		## For cn = 0, BAF~Unif(0,1)
		gamma0 <- dunif(bf, 0,1) * h[, 1]
		gammaTotal <- function(g) {
			sumg <- 0
			for(i in seq_len(length(g))) sumg <- sumg+g[[i]]
			sumg
		}
		##if(CHR=="chr10") browser()
		## BAF means for cn 1
		dA <- trNormal(mus["A"], sigmas["A"])
		dB <- trNormal(mus["B"], sigmas["B"])
		gamma1 <- calculateGamma(list(dA, dB), h[, 2], allele.prob[, c("A", "B")])
		m1 <- updateBMoments(gamma1, nms=c("A", "B"))

		## BAF means for cn 2
		dAB <- trNormal(mus["AB"], sigmas["AB"])
		gamma2 <- calculateGamma(list(dA, dAB, dB), h[, 3], allele.prob[, c("AA", "AB", "BB")])
		m2 <- updateBMoments(gamma2, nms=c("A", "AB", "B"))["AB", , drop=FALSE]

		## BAF means for cn 3
		dAAB <- trNormal(mus["AAB"], sigmas["AAB"])
		dABB <- trNormal(mus["ABB"], sigmas["ABB"])
		gamma3 <- calculateGamma(list(dA, dAAB, dABB, dB), h[, 5], allele.prob[, c("AAA", "AAB", "ABB", "BBB")])
		m3 <- updateBMoments(gamma3, nms=c("A", "AAB", "ABB", "B"))[c("AAB", "ABB"), ]

		dAAAB <- trNormal(mus["AAAB"], sigmas["AAAB"])
		dABBB <- trNormal(mus["ABBB"], sigmas["ABBB"])
		gamma4 <- calculateGamma(list(dA, dAAAB, dAB, dABBB, dB), h[, 6], allele.prob[,c("AAAA", "AAAB", "AABB", "ABBB", "BBBB")])
		m4 <- updateBMoments(gamma4, nms=c("A", "AAAB", "AB", "ABBB", "B"))[c("AAAB", "ABBB"), ]

		mus.new <- c(m1[,1], m2[, 1], m3[, 1], m4[, 1])
		mus.new["A"] <- 0
		mus.new["B"] <- 1
		names(mus.new)[3] <- "AB"
		if(anyMissing(mus.new))  {
			nas <- is.na(mus.new)
			mus.new[nas] <- theoreticalMeans[nas]
		}
		mus.new <- makeMusBafNondecreasing(mus.new)
		mus.new <- constrainMuBaf(mus.new)
		sigmas.new <- c(m1[,2], m2[, 2], m3[,2], m4[,2])
		names(sigmas.new)[3] <- "AB"
		sigmas.new <- sigmas.new[nms]
		sigmas.new <- pmax(sigmas.new, max(max(sigmas.new["A"], sigmas.new["B"]), 5e-4))

		## Update mixture components.
		Gtotal[, 1] <- gammaTotal(gamma1)
		Gtotal[, 2] <- gammaTotal(gamma2)
		Gtotal[, 3] <- gammaTotal(gamma3)
		Gtotal[, 4] <- gammaTotal(gamma4)
		## proportion outliers (how to account for homozygous deletions?)
		##    -- exclude consecutive outliers
		has.nas <- anyMissing(Gtotal)
		##nas <- is.na(Gtotal[,1])
		##has.nas <- any(nas)
		if(has.nas) {
			nas <- is.na(Gtotal[,1])
			Gtotal <- Gtotal[-which(nas), ]
		}
		isout <- rowMax(Gtotal) < 0.9
		p <- mean(diff(isout) != 0)
		params[["mus"]] <- mus.new
		params[["sigmas"]] <- sigmas.new
		##params[["sigmas"]] <- sigmas
		pOut[["prOutlier"]] <- p
		params[["prOutlier"]] <- pOut
		params
	}

	cn.sd.new <- rep(NA, S)
	updateCnParams <- function(cn, params, fv, bv, j) {
		mu <- params[["mu"]]
		sigma <- params[["sigmas"]]
		p <- params[["p"]]; q <- 1-p
		min.sd <- sds[j]/2

		h <- matrix(fv*bv, nr, S)
		m <- rowSums(h, na.rm=TRUE)
		h <- h/m

		d.unif <- cn.dunif[, j]
		d <- q*dnorm(cn, mean=matrix(mu, nr, S, byrow=TRUE), sd=matrix(sigma, nr, S, byrow=TRUE)) + p*d.unif
		d.total <- d+d.unif
		d1 <- d/(d+d.unif)
		d2 <- 1-d1
		##pout <- mean(rowMax(d1) < 0.5) ## uniform distribution has higher probability than any of the cn states

		g1 <- h * d1  ## emission from observed states
		g2 <- h * (1-d1) ## emission from outlier state

		total.g1 <- apply(g1, 2, sum, na.rm=TRUE)

		mu.new <- colSums(g1*cn, na.rm=TRUE)/total.g1
		mu.new[-4] <- constrainMu(mu.new[-4], is.log)
		mu.new[4] <- mu.new[normalIndex]
		mu.new <- makeNonDecreasing(mu.new)

		## For loop
		for(s in seq_len(S)) {
			cn.sd.new[s] <- sqrt(sum(g1[, s]*(cn-mu.new[s])^2, na.rm=TRUE)/total.g1[s])
		}
		cn.sd.new <- pmax(cn.sd.new, min.sd)
		cn.sd.new <- constrainSd(cn.sd.new)
		params[["mu"]] <- mu.new
		params[["sigmas"]] <- cn.sd.new

		total.g2 <- apply(g2, 2, sum, na.rm=TRUE)
		denom.eq52 <- apply(g1 + g2, 2, sum, na.rm=TRUE)
		p <- min(total.g2/denom.eq52)  ## altered states will have higher values of total.g2/denom.eq52.
		params[["p"]] <- min(p, 0.05)
		params
	}

	transitionPr <- function(TAUP, tauMAX=1-5e-6){
		p <- exp(-2*diff(position)/TAUP)
		minimum <- 1-1/((S-1)) + 0.01
		p <- pmax(p, minimum)
		p <- pmin(p, tauMAX)
		return(as.matrix(p))
	}
	tau <- transitionPr(TAUP)

	arm <- integer(nr)
	statePath <- integer(nr)
	scaleFactor <- rep(0, nr)
	bv <- fv <- matrix(0.0, nr*S, 1L)
	initialProb <- rep(1/S, S)

	fitViterbi <- function(emit){
		emit <- as.matrix(as.numeric(emit))
		res <- .C("viterbi2", emit=emit, pi=initialProb, tau=tau,
			  arm=arm, S=S, nr=nr, statePath=statePath,
			  fv=fv, bv=bv, 1, 1, 1, 3L, scaleFactor=scaleFactor)[c("statePath", "fv", "bv", "scaleFactor")]
		statePath <- res[["statePath"]]
		sf <- res[["scaleFactor"]]
		fv <- res[["fv"]]
		bv <- res[["bv"]]
		res <- list(statePath=statePath, fv=fv, bv=bv, sf=sf)
	}

	## computeLLR
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
	toGRanges <- function(statePath, j){
		id <- colnames(r)[j]
		rl <- Rle(statePath)
		starts <- position[start(rl)]
		ends <- position[end(rl)]
		states <- statePath[start(rl)]
		GRanges(CHR, IRanges(starts, ends), numberProbes=width(rl), state=states, sample=id)
	}

	if(computeLLR){
		log.initial <- log(initialProb)
		tauC <- 1-tau
		lP.N2N <- log(1-(tauC*(S-1))) ##probability normal -> normal
		lP.N2A <- log(tauC) ##probability normal -> altered
		P.A2A <- sapply(1-(tauC*(1+(S-2))), function(x) max(x, 0.01))
		lP.A2A <- log(P.A2A) ## probability altered to same altered state
		lP.A2N <- lP.N2A ##probability altered -> normal
		lP.A2Astar <- lP.N2A ## probability altered -> different altered state
		## featureRanges
		fr <- GRanges(rep(CHR, length(position)), IRanges(position, width=1))
	}
	computeLogLikRatio <- function(gr, emit){
		gr2 <- gr
		alt.index <- which(state(gr) != normalIndex & numberProbes(gr) > 1)
		if(length(alt.index)==0) return(rep(0.0, length(gr)))
		gr <- gr[alt.index, ]
		log.emission <- log(emit)
		L <- length(gr)
		LLR <- rep(NA,  L)
		olaps <- findOverlaps(gr, fr)
		index <- subjectHits(olaps)
		indexList <- split(index, queryHits(olaps))
		starts <- sapply(indexList, min)
		ends <- sapply(indexList, max)
		statePath <- as.integer(state(gr))
		T <- length(statePath)
		rangeLogLik <- function(from, to, thisState){
			index <- seq(from, to)
			index2 <- index[-1]## t=2, ...., t*
			if(from == 1){
				if(to < T){
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[to+1] + log.emission[to+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[to+1] + log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
				}
			} else { ## first index > 1
				index2 <- index[-1]
				if(to < T){
					logLik.vit <- lP.N2A[from] +
						sum(lP.A2A[index2]) +
							lP.A2N[to+1] +
								log.emission[from, thisState] +
									sum(log.emission[index2, thisState]) +
										log.emission[to+1, normalIndex]
					logLik.null <-
						sum(lP.N2N[index]) +
							lP.N2N[to+1]  +
								sum(log.emission[index, normalIndex]) +
									log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- lP.N2A[from] + log.emission[from, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
					logLik.null <- lP.N2N[from] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
				}
			}
			LLR <- logLik.vit - logLik.null
			return(LLR)
		}
		states <- state(gr)
		for(k in seq_len(L)) LLR[k] <- rangeLogLik(from=starts[k], to=ends[k], thisState=states[k])
		res <- rep(0.0, length(gr2))
		res[alt.index] <- LLR
		return(res)
	}

	list(initialCnParams=initialCnParams,
	     initialBafParams=initialBafParams,
	     updateBafEmission=updateBafEmission,
	     updateCnEmission=updateCnEmission,
	     updateBafParams=updateBafParams,
	     updateCnParams=updateCnParams,
	     trPr=transitionPr,
	     fitViterbi=fitViterbi,
	     toGRanges=toGRanges,
	     computeLogLikRatio=computeLogLikRatio)
}

getPrB <- function(pB=NULL){
	##n <- length(x)
	if(is.null(pB)){
		pB <- 0.5
	}
	pA <- 1-pB
	pAA <- pA^2
	pAB <- 2*pA*pB
	pBB <- 1-pAA-pAB
	pAAA <- pA^3
	pAAB <- 3*pA^2*pB
	pABB <- 3*pA*pB^2
	pBBB <- 1-pAAA-pAAB-pABB
	pAAAA <- pA^4
	pAAAB <- 4*pA^3*pB
	pAABB <- 6*pA^2*pB^2
	pABBB <- 4*pA*pB^3
	pBBBB <- pB^4
	p <- matrix(c(pB,
		      pA,
		      pAA,
		      pAB,
		      pBB,
		      pAAA,
		      pAAB,
		      pABB,
		      pBBB,
		      pAAAA,
		      pAAAB,
		      pAABB,
		      pABBB,
		      pBBBB),
		    byrow=FALSE, nrow=length(pB), ncol=14)
	colnames(p) <- c("B", "A", "AA", "AB", "BB",
			 "AAA", "AAB", "ABB", "BBB",
			 "AAAA", "AAAB", "AABB", "ABBB", "BBBB")
	return(p)
}

generatorFunG <- function(r, gt, is.snp, cnStates,
			  normalIndex, TAUP, limits, center,
			  position, is.log, computeLLR, chrom){
	S <- length(cnStates)
	nc <- ncol(r)
	nr <- nrow(r)
	np.index <- which(!is.snp)
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")

	## params for copy number
	sds <- apply(r, 2, mad, na.rm=TRUE)
	if(center) r <- centerCopyNumber(r, is.snp) + cnStates[normalIndex]
	r <- thresholdCopyNumber(r, limits=limits)

	initialCnParams <- function(j){
		mus <- cnStates
		sigmas <- rep(sds[j], S)
		##p <- matrix(c(0.99,  0.01), S, 2, byrow=TRUE)
		p <- 0.01
		paramsCN <- list(mu=mus, sigmas=sigmas, p=p)
	}

	g2 <- g1 <- matrix(NA, nr, S)
	d1 <- matrix(NA, nr, S)
	emitr <- matrix(NA, nr, S)
	cn.dunif <- dunif(r, limits[1], limits[2])
	Gtotal <- matrix(NA, nr, 4)

	## might try dt here instead
	mydnorm <- function(x){
		function(mean, sd){
			mus <- matrix(mean, nr, S, byrow=TRUE)
			sds <- matrix(sd, nr, S, byrow=TRUE)
			dnorm(x, mus, sds)
		}
	}

	updateCnEmission <- function(cn, params, j){
		mus <- params[["mu"]]
		sds <- params[["sigmas"]]
		p <- params[["p"]]
		q <- 1-p
		normalDens <- mydnorm(cn)
		d <- normalDens(mus, sds)
		emitr <- q*d + p*cn.dunif[, j]
		emitr.na <- is.na(emitr)
		hasna <- any(emitr.na)
		if(hasna) emitr[emitr.na] <- 1
		return(emitr)
	}
	cn.sd.new <- rep(NA, S)
	updateCnParams <- function(cn, params, fv, bv, j) {
		mu <- params[["mu"]]
		sigma <- params[["sigmas"]]
		p <- params[["p"]]; q <- 1-p
		min.sd <- sds[j]/2

		h <- matrix(fv*bv, nr, S)
		m <- rowSums(h, na.rm=TRUE)
		h <- h/m

		d.unif <- cn.dunif[, j]
		d <- q*dnorm(cn, mean=matrix(mu, nr, S, byrow=TRUE), sd=matrix(sigma, nr, S, byrow=TRUE)) + p*d.unif
		d1 <- d/(d+d.unif)

		g1 <- h * d1
		g2 <- h * (1-d1)

		total.g1 <- apply(g1, 2, sum, na.rm=TRUE)
		mu.new <- colSums(g1*cn, na.rm=TRUE)/total.g1
		mu.new[-4] <- constrainMu(mu.new[-4], is.log)
		mu.new[4] <- mu.new[normalIndex]
		mu.new <- makeNonDecreasing(mu.new)

		## For loop
		for(s in seq_len(S))  cn.sd.new[s] <- sqrt(sum(g1[, s]*(cn-mu.new[s])^2, na.rm=TRUE)/total.g1[s])
		cn.sd.new <- pmax(cn.sd.new, min.sd)
		cn.sd.new <- constrainSd(cn.sd.new)
		params[["mu"]] <- mu.new
		params[["sigmas"]] <- cn.sd.new

		## update p as the proportion of values not fit by any of the states
		total.g2 <- apply(g2, 2, sum, na.rm=TRUE)
		denom.eq52 <- apply(g1 + g2, 2, sum, na.rm=TRUE)
		p <- min(total.g2/denom.eq52)  ## altered
		params[["p"]] <- p
		params
	}

	transitionPr <- function(TAUP, tauMAX=1-5e-6){
		p <- exp(-2*diff(position)/TAUP)
		minimum <- 1-1/((S-1)) + 0.01
		p <- pmax(p, minimum)
		p <- pmin(p, tauMAX)
		return(as.matrix(p))
	}
	tau <- transitionPr(TAUP)

	arm <- integer(nr)
	statePath <- integer(nr)
	scaleFactor <- rep(0, nr)
	bv <- fv <- matrix(0.0, nr*S, 1L) ## forward and backward variables
	initialProb <- rep(1/S, S)

	fitViterbi <- function(emit){
		emit <- as.matrix(as.numeric(emit))
		res <- .C("viterbi2", emit=emit, pi=initialProb, tau=tau,
			  arm=arm, S=S, nr=nr, statePath=statePath,
			  fv=fv, bv=bv, 1, 1, 1, 3L,
			  scaleFactor=scaleFactor)[c("statePath", "fv", "bv", "scaleFactor")]
		statePath <- res[["statePath"]]
		sf <- res[["scaleFactor"]]
		fv <- res[["fv"]]
		bv <- res[["bv"]]
		res <- list(statePath=statePath, fv=fv, bv=bv, sf=sf)
	}

	## computeLLR
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
	toGRanges <- function(statePath, j){
		id <- colnames(r)[j]
		rl <- Rle(statePath)
		starts <- position[start(rl)]
		ends <- position[end(rl)]
		states <- statePath[start(rl)]
		GRanges(CHR, IRanges(starts, ends), numberProbes=width(rl), state=states, sample=id)
	}

	if(computeLLR){
		log.initial <- log(initialProb)
		tauC <- 1-tau
		lP.N2N <- log(1-(tauC*(S-1))) ##probability normal -> normal
		lP.N2A <- log(tauC) ##probability normal -> altered
		P.A2A <- sapply(1-(tauC*(1+(S-2))), function(x) max(x, 0.01))
		lP.A2A <- log(P.A2A) ## probability altered to same altered state
		lP.A2N <- lP.N2A ##probability altered -> normal
		lP.A2Astar <- lP.N2A ## probability altered -> different altered state
		## featureRanges
		fr <- GRanges(rep(CHR, length(position)), IRanges(position, width=1))
	}
	computeLogLikRatio <- function(gr, emit){
		gr2 <- gr
		alt.index <- which(state(gr) != normalIndex & numberProbes(gr) > 1)
		if(length(alt.index)==0) return(rep(0.0, length(gr)))
		gr <- gr[alt.index, ]
		log.emission <- log(emit)
		L <- length(gr)
		LLR <- rep(NA,  L)
		olaps <- findOverlaps(gr, fr)
		index <- subjectHits(olaps)
		indexList <- split(index, queryHits(olaps))
		starts <- sapply(indexList, min)
		ends <- sapply(indexList, max)
		statePath <- as.integer(state(gr))
		T <- length(statePath)
		rangeLogLik <- function(from, to, thisState){
			index <- seq(from, to)
			index2 <- index[-1]## t=2, ...., t*
			if(from == 1){
				if(to < T){
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[to+1] + log.emission[to+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[to+1] + log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
				}
			} else { ## first index > 1
				index2 <- index[-1]
				if(to < T){
					logLik.vit <- lP.N2A[from] +
						sum(lP.A2A[index2]) +
							lP.A2N[to+1] +
								log.emission[from, thisState] +
									sum(log.emission[index2, thisState]) +
										log.emission[to+1, normalIndex]
					logLik.null <-
						sum(lP.N2N[index]) +
							lP.N2N[to+1]  +
								sum(log.emission[index, normalIndex]) +
									log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- lP.N2A[from] + log.emission[from, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
					logLik.null <- lP.N2N[from] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
				}
			}
			LLR <- logLik.vit - logLik.null
			return(LLR)
		}
		states <- state(gr)
		for(k in seq_len(L)) LLR[k] <- rangeLogLik(from=starts[k], to=ends[k], thisState=states[k])
		res <- rep(0.0, length(gr2))
		res[alt.index] <- LLR
		return(res)
	}
	list(initialCnParams=initialCnParams,
	     updateCnEmission=updateCnEmission,
	     updateCnParams=updateCnParams,
	     trPr=transitionPr,
	     fitViterbi=fitViterbi,
	     toGRanges=toGRanges,
	     computeLogLikRatio=computeLogLikRatio)
}



gcPercent <- function(feature.gr, width){
	##feature.gr <- feature.gr[chromosome(feature.gr) == chrom, ]
	bsgenomepkg <- paste("BSgenome.Hsapiens.UCSC.", metadata(feature.gr)$genome, sep="")
	library(bsgenomepkg, character.only=TRUE)
	Hsapiens <- get("Hsapiens")
	chrom <- chromosome(feature.gr)[1]
	chr <- Hsapiens[[chrom]]
	v <- suppressWarnings(Views(chr, start=start(feature.gr), end=end(feature.gr)))
	gcFunction <- function(x){
		alf <- alphabetFrequency(x, as.prob=TRUE)
		rowSums(alf[, c("G", "C")])
	}
	probeGC <- suppressWarnings(gcFunction(v))
	probeGC <- as.integer(probeGC*100)
	return(probeGC)
}

gcCorrect <- function(object, width){
	feature.gr <- makeFeatureGRanges(object)
	feature.gr <- sort(feature.gr)
	feature.gr <- resize(feature.gr, width=width, fix="center")

	arm <- .getArm(chromosome(object), position(object), genomeBuild(object))
	arm <- factor(arm, levels=unique(arm))
	feature.grl <- split(feature.gr, arm)

	gcP <- lapply(feature.grl, gcPercent, width=width)
	binGCcontent <- function(x) bins <- cut(x, breaks=seq(0,100), labels=FALSE)
	binList <- lapply(gcP, binGCcontent)

	## Calculate mean lrr for each bin.
	## Subtract the mean lrr for the gc bins.  Add back the grand median so that the
	## center for the chromosome is unchanged
	avgRBin <- function(r, bin, grand.median) {
		bins <- as.integer(factor(bin))
		ravg <- sapply(split(r, bins), mean, na.rm=TRUE)
		## what about QN by middle bins..
		ravg <- ravg[bins]
		r.adj <- r - ravg + grand.median
		##boxplot(split(r/100, bins), pch=".")
		##boxplot(split(r.adj/100, bins), pch=".")
		return(r.adj)
	}
	for(j in seq_len(ncol(object))){
		rList <- split(lrr(object)[, j], arm)  ## can do this on the integer scale
		medianR <- sapply(rList, median, na.rm=TRUE)
		adjR <- mapply(avgRBin, rList, binList, medianR)
		rvec <- unlist(adjR)
		lrr(object)[, j] <- rvec
	}
	return(object)
}


