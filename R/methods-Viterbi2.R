##setMethod("initialize", signature(.Object="Viterbi2"),
##	  function(.Object,
##		   numberFeatures=0L,
##		   numberStates=6L,
##		   emission=matrix(0L, numberFeatures*numberStates, 1L),
##		   initialProb=rep(1/numberStates, numberStates),
##		   arm=integer(numberFeatures),
##		   transitionProb=matrix(numeric(max(c(numberFeatures-1, 0))), ncol=1L),
##		   viterbiStatePath=integer(numberFeatures),
##		   normalIndex=3L,
##		   forwardVariable=matrix(0, nrow(emission), 1L),
##		   backwardVariable=matrix(0, nrow(emission), 1L),
##		   normal2altered=1,
##		   altered2normal=1,
##		   altered2altered=1,
##		   scaleFactor=rep(0, numberFeatures), ...){
##		  callNextMethod(.Object,
##				 numberFeatures=numberFeatures,
##				 numberStates=numberStates,
##				 emission=emission,
##				 initialProb=initialProb,
##				 arm=arm,
##				 transitionProb=transitionProb,
##				 viterbiStatePath=viterbiStatePath,
##				 normalIndex=normalIndex,
##				 forwardVariable=forwardVariable,
##				 backwardVariable=backwardVariable,
##				 normal2altered=normal2altered,
##				 altered2normal=altered2normal,
##				 altered2altered=altered2altered,
##				 scaleFactor=scaleFactor, ...)
##	  })
##
##setValidity("Viterbi2", function(object){
##	ns <- nStates(object)
##	##nf <- nrow(object)
##	nf <- object@numberFeatures
##	nf.ns <- ns*nf
##	emit <- emission(object)
##	if(ncol(emit) != 1){
##		return("emission matrix should have a single column")
##	}
##	if(nrow(emit) != nf.ns){
##		return("emission matrix should have number of features x number of states elements")
##	}
##	bv <- backwardVariable(object)
##	fv <- forwardVariable(object)
##	if(!all(dim(emit) == dim(bv))){
##		return("emission matrix and backwardVariable should have the same dimension")
##	}
##	if(!all(dim(emit) == dim(fv))){
##		return("emission matrix and forwardVariable should have the same dimension")
##	}
##	statePath <- viterbiStatePath(object)
##	if(length(statePath) != nf){
##		return("viterbiStatePath should have the same length as the number of features")
##	}
##	if(ns > 0){
##		ni <- normalStateIndex(object)
##		if(ni <= 0 | ni > ns){
##			return("the index for the normal state (normalIndex) \n must be an integer greater than zero and \n less than or equal to the number of states")
##		}
##		if(sum(initialProb(object)) != 1)
##			return("initialProb must sum to 1")
##	}
##	if(nf > 0){
##		tp <- transitionProb(object)
##		if(length(tp) != (nf-1)){
##			return("transitionProb vector should be the number of features -1")
##		}
##	}
##})


##setMethod("backwardVariable", signature(object="Viterbi2"),
##	  function(object){
##		  object@backwardVariable
##	  })
##setMethod("forwardVariable", signature(object="Viterbi2"),
##	  function(object){
##		  object@forwardVariable
##	  })
##
##
##setMethod("scaleFactor", signature(object="Viterbi2"),
##	  function(object){
##		  object@scaleFactor
##	  })
##
##setMethod("fit", signature(object="Viterbi2"), function(object){
##	fitViterbi2(object)
##})
##
##fitViterbi2 <- function(object){
##	tmp <- .C("viterbi2",
##		  emission(object),
##		  initialProb(object),
##		  transitionProb(object),
##		  object@arm,
##		  nStates(object),
##		  object@numberFeatures,
##		  viterbiStatePath(object),
##		  forwardVariable(object),
##		  backwardVariable(object),
##		  object@normal2altered,
##		  object@altered2normal,
##		  object@altered2altered,
##		  normalStateIndex(object),
##		  scaleFactor(object))
##	new("Viterbi2",
##	    emission=emission(object),
##	    initialProb=initialProb(object),
##	    transitionProb=transitionProb(object),
##	    arm=object@arm,
##	    numberStates=nStates(object),
##	    numberFeatures=object@numberFeatures,
##	    viterbiStatePath=tmp[[7]],
##	    forwardVariable=tmp[[8]],
##	    backwardVariable=tmp[[9]],
##	    normal2altered=object@normal2altered,
##	    altered2normal=object@altered2normal,
##	    altered2altered=object@altered2altered,
##	    normalIndex=normalStateIndex(object),
##	    scaleFactor=tmp[[14]])
##}
##
##
##setMethod("show", signature(object="Viterbi2"), function(object){
##	cat("Number of features:", object@numberFeatures, "\n")
##	cat("Number of states:", nStates(object), "\n")
##	loglik <- -round(sum(log(scaleFactor(object))))
##	if(is.finite(loglik))
##		cat("loglik:", loglik, "\n")
##})






##		for(j in seq_len(nc)){
##			tmp <- new("Viterbi2",
##				   emission=as.matrix(as.numeric(e[, j, ])),
##				   initialProb=initialProb,
##				   transitionProb=taus,
##				   numberStates=length(cnStates),
##				   numberFeatures=nr,
##				   normalIndex=normalIndex)
##			viterbiList[[j]] <- fitViterbi2(tmp)
##			loglik[i,j] <- -sum(log(scaleFactor(viterbiList[[j]])))
##		}
		## fit viterbi
##	paramsBaf <- updateFun$bafUp(paramsBaf, fv, bv)
##	paramsCn <- updateFun$cnUp(paramsCn, fv, bv)
##	emitB <- updateFun$bafE()
##	emitCn <- updateFun$cnE(paramsCn)
#	}
##	for(i in seq_len(nupdates)){
##		mu.sigma <- list()
##		for(j in seq_len(J)) mu.sigma[[j]] <- updateMu(x=r[, j], mu=mus[j, ], sigma=sigmas[j, ], normalIndex=normalIndex, rohIndex=rohIndex, limits=limits,  fv=forwardVariable(viterbiList[[j]]),  bv=backwardVariable(viterbiList[[j]]), constrainMu=constraint[["mus.cn"]], is.log=is.log)
##		mus <- lapply(mu.sigma, "[[", 1)
##		mus <- do.call("rbind", mus)
##		sigmas <- lapply(mu.sigma, "[[", 2)
##		sigmas <- do.call("rbind", sigmas)
##		p <- lapply(mu.sigma, "[[", 3)
##		if(!missing(b)){
##			outlierProb <- prOutlierBAF[["initial"]]
##			mu.sigmaBAF <- list()
##			for(j in seq_len(J)){
##				prOutlierBAF[["initial"]] <- outlierProb[j, ]
##				mu.sigmaBAF[[j]] <- updateBaf(x=b[, j], cnStates=cnStates, mus=musBAF[j, ], sigmas=sdsBAF[j, ], normalIndex=normalIndex,  fv=forwardVariable(viterbiList[[j]]), bv=backwardVariable(viterbiList[[j]]), prOutlierBAF=prOutlierBAF, nUpdates=1, constraint=constraint[2:3], allele.prob=allele.prob)
##			}
##			musBAF <- do.call("rbind", lapply(mu.sigmaBAF, "[[", 1))
##			sdsBAF <- do.call("rbind", lapply(mu.sigmaBAF, "[[", 2))
##			## update for next iteration
##			prOutlierBAF[["initial"]] <- do.call("rbind", lapply(mu.sigmaBAF, "[[", 3)) ## outlier probability
##		}
##		emitr <- cnEmissionFromMatrix(object=r, mus=mus, sds=sigmas, p=p, limits=limits, log.it=FALSE, ...)
##		if(!missing(b)) emitb <- bafEmissionFromMatrix2(BAF=b, is.snp=is.snp, mus=musBAF, sds=sdsBAF, prOutlier=prOutlierBAF[["initial"]], p.hom=p.hom, log.it=FALSE, allele.prob=allele.prob,...)
##		e <- emitb*emitr
##		viterbiList <- list()
##		for(j in seq_len(J)){
##			tmp <- new("Viterbi2", emission=as.matrix(as.numeric(e[, j, ])), initialProb=initialProb, transitionProb=taus, numberStates=length(cnStates), numberFeatures=nr, normalIndex=normalIndex)
##			viterbiList[[j]] <- fitViterbi2(tmp)
##			loglik[i,j] <- -sum(log(scaleFactor(viterbiList[[j]])))
##		}
##		if(i > 1){
##			llr <- loglik[i, ]-loglik[i-1, ]
##			check.break <- all(abs(llr) < tolerance)
##		} else check.break <- FALSE
##		if(check.break) break()
##	}
##	loglik <- loglik[seq_len(i), ]
##	if(!returnViterbiObject){
##		id <- object <- NULL
##		res <- foreach(object=viterbiList, id=colnames(r)) %do% {
##			computeLoglikFromViterbi2(object=object,
##						  chrom=chrom,
##						  id=id,
##						  pos=pos)
##		}
##		res <- GRangesList(res)
##		##res <- stackRangedData(res)
##		##res <- GRangesList(res)
##		names(res) <- colnames(r)
##		return(res)
##		##return(list(rangedData=res, loglik=loglik))
##	} else {
##		return(viterbiList)
##	}
##}
