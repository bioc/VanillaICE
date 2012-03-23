setMethod("initialize", signature(.Object="Viterbi2"),
	  function(.Object,
		   numberFeatures=0L,
		   numberStates=6L,
		   emission=matrix(0L, numberFeatures*numberStates, 1L),
		   initialProb=rep(1/numberStates, numberStates),
		   arm=integer(numberFeatures),
		   transitionProb=matrix(numeric(max(c(numberFeatures-1, 0))), ncol=1L),
		   viterbiStatePath=integer(numberFeatures),
		   normalIndex=3L,
		   forwardVariable=matrix(0, nrow(emission), 1L),
		   backwardVariable=matrix(0, nrow(emission), 1L),
		   normal2altered=1,
		   altered2normal=1,
		   altered2altered=1,
		   scaleFactor=rep(0, numberFeatures), ...){
		  callNextMethod(.Object,
				 numberFeatures=numberFeatures,
				 numberStates=numberStates,
				 emission=emission,
				 initialProb=initialProb,
				 arm=arm,
				 transitionProb=transitionProb,
				 viterbiStatePath=viterbiStatePath,
				 normalIndex=normalIndex,
				 forwardVariable=forwardVariable,
				 backwardVariable=backwardVariable,
				 normal2altered=normal2altered,
				 altered2normal=altered2normal,
				 altered2altered=altered2altered,
				 scaleFactor=scaleFactor, ...)
	  })

setValidity("Viterbi2", function(object){
	ns <- nStates(object)
	##nf <- nrow(object)
	nf <- object@numberFeatures
	nf.ns <- ns*nf
	emit <- emission(object)
	if(ncol(emit) != 1){
		return("emission matrix should have a single column")
	}
	if(nrow(emit) != nf.ns){
		return("emission matrix should have number of features x number of states elements")
	}
	bv <- backwardVariable(object)
	fv <- forwardVariable(object)
	if(!all(dim(emit) == dim(bv))){
		return("emission matrix and backwardVariable should have the same dimension")
	}
	if(!all(dim(emit) == dim(fv))){
		return("emission matrix and forwardVariable should have the same dimension")
	}
	statePath <- viterbiStatePath(object)
	if(length(statePath) != nf){
		return("viterbiStatePath should have the same length as the number of features")
	}
	if(ns > 0){
		ni <- normalStateIndex(object)
		if(ni <= 0 | ni > ns){
			return("the index for the normal state (normalIndex) \n must be an integer greater than zero and \n less than or equal to the number of states")
		}
		if(sum(initialProb(object)) != 1)
			return("initialProb must sum to 1")
	}
	if(nf > 0){
		tp <- transitionProb(object)
		if(length(tp) != (nf-1)){
			return("transitionProb vector should be the number of features -1")
		}
	}
})


setMethod("backwardVariable", signature(object="Viterbi2"),
	  function(object){
		  object@backwardVariable
	  })
setMethod("forwardVariable", signature(object="Viterbi2"),
	  function(object){
		  object@forwardVariable
	  })


setMethod("scaleFactor", signature(object="Viterbi2"),
	  function(object){
		  object@scaleFactor
	  })

setMethod("fit", signature(object="Viterbi2"), function(object){
	fitViterbi2(object)
})

fitViterbi2 <- function(object){
	tmp <- .C("viterbi2",
		  emission(object),
		  initialProb(object),
		  transitionProb(object),
		  object@arm,
		  nStates(object),
		  object@numberFeatures,
		  viterbiStatePath(object),
		  forwardVariable(object),
		  backwardVariable(object),
		  object@normal2altered,
		  object@altered2normal,
		  object@altered2altered,
		  normalStateIndex(object),
		  scaleFactor(object))
	new("Viterbi2",
	    emission=emission(object),
	    initialProb=initialProb(object),
	    transitionProb=transitionProb(object),
	    arm=object@arm,
	    numberStates=nStates(object),
	    numberFeatures=object@numberFeatures,
	    viterbiStatePath=tmp[[7]],
	    forwardVariable=tmp[[8]],
	    backwardVariable=tmp[[9]],
	    normal2altered=object@normal2altered,
	    altered2normal=object@altered2normal,
	    altered2altered=object@altered2altered,
	    normalIndex=normalStateIndex(object),
	    scaleFactor=tmp[[14]])
}

##.checksViterbi <- function(initialCnStateMeans,
##			   normalIndex,
##			   is.log){
##	if(any(diff(initialCnStateMeans) < 0)) stop("initialCnStateMeans must be ordered")
##	if(missing(is.log)) stop("is.log can not be missing. \n Must specify whether the copy number estimates are on the log-scale (T/F). ")
##	if(missing(normalIndex)) stop("normalIndex can not be missing")
##	if(!normalIndex %in% seq_along(initialCnStateMeans)){
##		stop("normalIndex must be greater than zero and less than or equal to the number of states")
##	}
##}


##mmFromOligoSnpSet <- function(object,
##			       initialCnStateMeans=c(0, 1, 2, 2, 3, 4),
##			       normalIndex,
##			       TAUP=1e8,
##			       is.log,
##			       initialProb=rep(1/length(initialCnStateMeans), length(initialCnStateMeans)),
##			       ...){
##	.checksViterbi(initialCnStateMeans=initialCnStateMeans,
##		       normalIndex=normalIndex,
##		       is.log=is.log)
##	S <- length(cnStates)
##	object <- chromosomePositionOrder(object)
##	arm <- .getArm(chromosome(object), position(object))
##	indices <- split(seq_along(arm), arm)
##
##	## for each arm
##	##      - compute initial emission probability estimates
##	##      - initialize Viterbi object
##	##      - run viterbi
##	##      - compute log lik
##	##      - update model parameters and emission probs.
##	##      - update emission probs in viterbi object
##	##      - run viterbi
##	##      - compare log lik to previous log lik
##	##        ...
##	chrlist <- split(chromosome(object), arm)
##	poslist <- split(position(object), arm)
##	armlist <- split(arm, arm)
##	snplist <- split(isSnp(object), arm)
##
##	taulist <- foreach(x=poslist) %do% computeTransitionProb(x=x,
##			   TAUP=TAUP, S=S)
##	##taulist <- lapply(taulist, as.matrix)
##	nf <- sapply(armlist, length)
##	if(use.baf){
##		baflist <- split(baf(object), arm)
##		baflist <- lapply(baflist, as.matrix)
##	} else {
##		gtlist <- split(calls(object), arm)
##		gtlist <- lapply(baflist, as.matrix)
##	}
##	if(use.baf){
##		emitb <- foreach(object=baflist,
##				 is.snp=snplist,
##				 .packages="VanillaICE") %do% {
##					 bafEmission(object=object,
##						     is.snp=is.snp, ...)
##				 }
##	} else {
##		stop("emission for genotypes")
##	}
##	emitr <- foreach(object=lrrlist,
##			 is.snp=snplist,
##			 chrom=chrlist,
##			 .packages="VanillaICE") %do% {
##				 cnEmission(object=object,
##					    is.snp=is.snp,
##					    cnStates=initialCnStateMeans,
##					    is.log=is.log,
##					    normalIndex=3,...)
##
##			 }
##	emitlist <- foreach(r=emitr, b=emitb) %do% as.matrix(as.numeric(r[,1,]+b[,1,]))
##	objectlist <- foreach(arm=armlist,
##			      tau=taulist,
##			      emission=emitList,
##			      ) %do% {
##				      new("Viterbi2",
##					  numberFeatures=length(arm),
##					  numberStates=S,
##					  transitionProb=tau,
##					  initialProb=initialProb,
##					  normalIndex=normalIndex,
##					  emission=emission)
##			      }
##
##
##


setMethod("show", signature(object="Viterbi2"), function(object){
	cat("Number of features:", object@numberFeatures, "\n")
	cat("Number of states:", nStates(object), "\n")
	loglik <- -round(sum(log(scaleFactor(object))))
	if(is.finite(loglik))
		cat("loglik:", loglik, "\n")
})


viterbi2Wrapper <- function(r, b, gt, pos, is.snp, cnStates,
			    chrom,
			    prOutlierBAF=1e-3, p.hom=0.05, TAUP=1e8,
			    is.log, center=TRUE,
			    reestimation=TRUE,
			    limits,
			    initialProb=rep(1/length(cnStates), length(cnStates)),
			    normalIndex=3L,
			    rohIndex=normalIndex+1L,
			    nupdates=10,
			    tolerance=1,
			    returnViterbiObject=FALSE, ...){
	S <- length(cnStates)
	J <- ncol(r)
	## todo: reestimation of mixture probs for BAFs
	if(!missing(b)){
		emitb <- array(NA, dim=c(nrow(b), J, S))
		for(j in seq_len(J)){
			emitb[, j, ] <- bafEmission(object=b[, j, drop=FALSE],
						    is.snp=is.snp,
						    prOutlier=prOutlierBAF,
						    p.hom=p.hom, log.it=FALSE,...)[, 1, ]
		}
	} else {
		stopifnot(rohIndex==4L)
		emitb <- gtEmission(object=gt,
				    is.snp=is.snp,
				    rohIndex=c(2L, 4L),
				    prGtHom=c(2/3, 0.99, 0.7, 0.99, 1/2, 2/5),
				    S=S, log.it=FALSE)
	}
	taus <- computeTransitionProb(x=pos, TAUP=TAUP, S=S)
	sds <- apply(r, 2, mad, na.rm=TRUE)
	if(center){
		r <- centerCopyNumber(r, is.snp) + cnStates[normalIndex]
	}
	r <- thresholdCopyNumber(r, limits=limits)
	loglik <- matrix(NA, nupdates, J)
	viterbiList <- vector("list", J)
	for(i in seq_len(nupdates)){
		if(i == 1){
			mus <- matrix(cnStates, J, S, byrow=TRUE)
			sigmas <- matrix(sds, J, S, byrow=FALSE)
			## homozygous deletions often have a higher variance.
			## Initialize the variance for this state to have a bigger variance
			sigmas[, 1] <- sigmas[, 1] * 4
			sigmas[, 2] <- sigmas[, 2] * 2
			p <- replicate(J, matrix(c(0.99, 0.01), S, 2, byrow=TRUE), simplify=FALSE)
		} else {
			j <- NULL
			mu.sigma <- foreach(j = seq_len(J)) %do% {
				updateMu(x=r[, j],
					 mu=mus[j, ],
					 sigma=sigmas[j, ],
					 normalIndex=normalIndex,
					 rohIndex=rohIndex,
					 limits=limits,
					 is.log=is.log,
					 fv=forwardVariable(viterbiList[[j]]),
					 bv=backwardVariable(viterbiList[[j]]),
					 nUpdates=1)
			}
			mus <- lapply(mu.sigma, "[[", 1)
			mus <- do.call("rbind", mus)
			sigmas <- lapply(mu.sigma, "[[", 2)
			sigmas <- do.call("rbind", sigmas)
			p <- lapply(mu.sigma, "[[", 3)
		}
		emitr <- cnEmissionFromMatrix(object=r,
					      mus=mus,
					      sds=sigmas,
					      p=p,
					      limits=limits,
					      log.it=FALSE, ...)
		e <- emitb*emitr
		viterbiList <- list()
		for(j in seq_len(J)){
			tmp <- new("Viterbi2",
				   emission=as.matrix(as.numeric(e[, j, ])),
				   initialProb=initialProb,
				   transitionProb=taus,
				   numberStates=length(cnStates),
				   numberFeatures=nrow(r),
				   normalIndex=normalIndex)
			viterbiList[[j]] <- fitViterbi2(tmp)
			loglik[i,j] <- -sum(log(scaleFactor(viterbiList[[j]])))
		}
		if(i > 1){
			llr <- loglik[i, ]-loglik[i-1, ]
			if(all(llr < tolerance)){
				loglik <- loglik[1:i, ]
				break()
			}
		}
	}
	if(!returnViterbiObject){
		id <- object <- NULL
		res <- foreach(object=viterbiList, id=colnames(r)) %do% {
			computeLoglikFromViterbi2(object=object,
						  chrom=chrom, id=id,
						  pos=pos)
		}
		res <- stackRangedData(res)
		##return(list(viterbi=viterbiObject, loglik=loglik))
		return(list(rangedData=res, loglik=loglik))
	} else {
		return(viterbiList)
	}
}
