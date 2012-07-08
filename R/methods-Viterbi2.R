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


setMethod("show", signature(object="Viterbi2"), function(object){
	cat("Number of features:", object@numberFeatures, "\n")
	cat("Number of states:", nStates(object), "\n")
	loglik <- -round(sum(log(scaleFactor(object))))
	if(is.finite(loglik))
		cat("loglik:", loglik, "\n")
})


viterbi2Wrapper <- function(r, b, pos, is.snp, cnStates,
			    chrom,
			    prOutlierBAF=list(initial=1e-5, max=1e-3, maxROH=1e-5),
			    p.hom=0.05,
			    TAUP=1e8,
			    is.log,
			    center=TRUE,
			    limits,
			    initialProb=rep(1/length(cnStates), length(cnStates)),
			    normalIndex=3L,
			    nupdates=10,
			    tolerance=5,
			    computeLLR=TRUE,
			    returnEmission=FALSE,
			    ...){
	nc <- ncol(r)
	updateFun <- generatorFun(r, b, is.snp=is.snp, cnStates=cnStates,
				  normalIndex=normalIndex,  TAUP=TAUP,
				  limits=limits, center=center,
				  prOutlierBAF=prOutlierBAF, p.hom=p.hom,
				  position=pos, is.log=is.log,
				  computeLLR=computeLLR, chrom=chrom)
	paramsB <- updateFun$initialBafParams()
	statePath <- matrix(NA, nrow(r), nc)
	grl <- vector("list", nc)
	if(returnEmission) emitArray <- array(NA, dim=c(nrow(r), nc, length(cnStates)))
	## Things that change with iterator i:
	##  --  first two moments
	##  --  mixture probabilities
	## pseudocode
	for(j in seq_len(nc)){
		paramsC <- updateFun$initialCnParams(j)
		## if we subset r and bf here we wouldn't need to do "[" for each update i
		bf <- b[, j]
		cn <- r[, j]
		sameStatePath <- 0
		llr <- loglik <- rep(0, nupdates)
		for(i in seq_len(nupdates)){
			paramsBprev <- paramsB
			paramsCprev <- paramsC
			emitB <- updateFun$updateBafEmission(bf, paramsB, j)
			emitC <- updateFun$updateCnEmission(cn, paramsC, j)
			emit <- emitB*emitC
			vit.res <- updateFun$fitViterbi(emit)
			if(i > 1) sameStatePath <- if(identical(statePath, vit.res[["statePath"]])) sameStatePath+1 else 0
			statePath <- vit.res[["statePath"]]
			fv <- vit.res[["fv"]]
			bv <- vit.res[["bv"]]
			paramsB <- updateFun$updateBafParams(bf, paramsB, fv, bv, j)
			paramsC <- updateFun$updateCnParams(cn, paramsC, fv, bv, j)
			## can we compare the scale factor across runs?
			## Equation 103, Rabiner
			loglik[i] <- sum(log(vit.res[["sf"]]))
			if(i == 1) next()
			llr[i] <- loglik[i]-loglik[i-1]
			cond <- (abs(llr[i]) < tolerance) || (sameStatePath > 1)
			if(cond | i == nupdates) {
				statePath <- vit.res[["statePath"]]
				gr <- updateFun$toGRanges(statePath, j)
				if(computeLLR) elementMetadata(gr)$LLR <- updateFun$computeLogLikRatio(gr, emit)
				grl[[j]] <- gr
				break()
			}
		}
		if(returnEmission) emitArray[, j, ] <- emit
	}
	if(returnEmission) return(emitArray)
	grl <- GRangesList(grl)
	names(grl) <- colnames(r)
	return(grl)
}

viterbiWrapperG <- function(r, gt, pos, is.snp, cnStates,
			     chrom,
			     p.hom=0.05,
			     TAUP=1e8,
			     is.log,
			     center=TRUE,
			     limits,
			     initialProb=rep(1/length(cnStates), length(cnStates)),
			     normalIndex=3L,
			     nupdates=10,
			     tolerance=5,
			     computeLLR=TRUE,
			     ...){
	emitG <- gtEmission(object=gt, is.snp=is.snp, rohIndex=c(2L, 4L), prGtHom=c(2/3, 0.99, 0.7, 0.99, 1/2, 2/5), S=length(cnStates), log.it=FALSE)
	nc <- ncol(r)
	updateFun <- generatorFunG(r, gt, is.snp=is.snp,
				   cnStates=cnStates, normalIndex=normalIndex,
				   TAUP=TAUP, limits=limits, center=center,
				   position=pos, is.log=is.log,
				   computeLLR=computeLLR, chrom=chrom)
	statePath <- matrix(NA, nrow(r), nc)
	grl <- vector("list", nc)
	## Things that change with iterator i:
	##  --  first two moments
	##  --  mixture probabilities
	## pseudocode
	for(j in seq_len(nc)){
		paramsC <- updateFun$initialCnParams(j)
		##gt <- g[, j]
		cn <- r[, j]
		sameStatePath <- 0
		llr <- loglik <- rep(0, nupdates)
		for(i in seq_len(nupdates)){
			paramsCprev <- paramsC
			emitC <- updateFun$updateCnEmission(cn, paramsC, j)
			emit <- emitG[,j,]*emitC
			vit.res <- updateFun$fitViterbi(emit)
			if(i > 1) sameStatePath <- if(identical(statePath, vit.res[["statePath"]])) sameStatePath+1 else 0
			statePath <- vit.res[["statePath"]]
			fv <- vit.res[["fv"]]
			bv <- vit.res[["bv"]]
			##paramsB <- updateFun$updateBafParams(bf, paramsB, fv, bv, j)
			paramsC <- updateFun$updateCnParams(cn, paramsC, fv, bv, j)
			loglik[i] <- sum(log(vit.res[["sf"]]))
			if(i == 1) next()
			llr[i] <- loglik[i]-loglik[i-1]
			cond <- abs(llr[i]) < tolerance || sameStatePath>1
			if(cond | i == nupdates) {
				statePath <- vit.res[["statePath"]]
				gr <- updateFun$toGRanges(statePath, j)
				if(computeLLR) elementMetadata(gr)$LLR <- updateFun$computeLogLikRatio(gr, emit)
				grl[[j]] <- gr
				break()
			}
		}
	}
	grl <- GRangesList(grl)
	names(grl) <- colnames(r)
	return(grl)
}

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
