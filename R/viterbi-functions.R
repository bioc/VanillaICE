viterbi2Wrapper <- function(index.samples,
			    cnStates,
			    prOutlierBAF=list(initial=1e-5, max=1e-3, maxROH=1e-5),
			    p.hom=0.05,
			    is.log,
			    limits,
			    normalIndex=3L,
			    nupdates=10,
			    tolerance=5,
			    computeLLR=TRUE,
			    returnEmission=FALSE,
			    verbose=FALSE,
			    grFun,
			    matrixFun,
			    snp.index,
			    anyNP){
	toGRanges <- grFun$toGRanges
	transitionProbs <- grFun$transitionProbs
	tau <- grFun$tau
	b.r <- matrixFun(index.samples)
	b <- b.r[["b"]]
	r <- b.r[["r"]]
	rm(b.r, grFun); gc()
	nc <- ncol(r); nr <- nrow(r); S <- length(cnStates)
	ids <- colnames(r)
	##if(center) r <- centerCopyNumber(r, is.snp) + cnStates[normalIndex]
	r <- thresholdCopyNumber(r, limits=limits)
	updateFun <- generatorFun(r, b,
				  snp.index=snp.index,
				  cnStates=cnStates,
				  normalIndex=normalIndex,
				  tau=tau,
				  limits=limits,
				  prOutlierBAF=prOutlierBAF,
				  p.hom=p.hom,
				  is.log=is.log,
				  computeLLR=computeLLR,
				  verbose=verbose,
				  transitionProbs=transitionProbs)
	statePath <- matrix(NA, nrow(r), nc)
	grl <- vector("list", nc)
	if(returnEmission) emitArray <- array(NA, dim=c(nrow(r), nc, length(cnStates)))
	## Things that change with iterator i:
	##  --  first two moments
	##  --  mixture probabilities
	if(anyNP) emitB <- matrix(1, nrow(r), length(cnStates))
	for(j in seq_len(nc)){
		paramsB <- updateFun$initialBafParams()
		paramsC <- updateFun$initialCnParams(j)
		## if we subset r and bf here we wouldn't need to do "[" for each update i
		##if(anyNP) {
		##bf <- b[snp.index, j]
		##} else bf <- b[, j]
		bf <- b[, j]
		cn <- r[, j]
		sameStatePath <- 0
		llr <- loglik <- rep(0, nupdates)
		for(i in seq_len(nupdates)){
			paramsBprev <- paramsB
			paramsCprev <- paramsC
			if(anyNP){
				emitB[snp.index, ] <- updateFun$updateBafEmission(bf, paramsB, j)
			} else emitB <- updateFun$updateBafEmission(bf, paramsB, j)
			emitC <- updateFun$updateCnEmission(cn, paramsC, j)
			emit <- emitB*emitC
			rm(emitC)
			vit.res <- updateFun$fitViterbi(emit)
			if(i > 1) sameStatePath <- if(identical(statePath, vit.res[["statePath"]])) sameStatePath+1 else 0
			statePath <- vit.res[["statePath"]]
			h <- matrix(vit.res[["fv"]]*vit.res[["bv"]],
				    nr, S)
			if(anyNP){
				paramsB <- updateFun$updateBafParams(bf, paramsB, h[snp.index, ], j)
			} else paramsB <- updateFun$updateBafParams(bf, paramsB, h, j)
			if(verbose) print(paramsB)
			paramsC <- updateFun$updateCnParams(cn, paramsC, h, j)
			if(verbose) print(paramsC)
			## can we compare the scale factor across runs?
			## Equation 103, Rabiner
			loglik[i] <- sum(log(vit.res[["sf"]]))
			if(verbose) message(loglik[i])
			##if(i == 1) next()
			if(i > 1){
				llr[i] <- loglik[i]-loglik[i-1]
				cond <- (abs(llr[i]) < tolerance) || (sameStatePath > 1)
			} else cond <- FALSE
			if(cond | i == nupdates) {
				if(!returnEmission){
					statePath <- vit.res[["statePath"]]
					##gr <- updateFun$toGRanges(statePath, j)
					gr <- toGRanges(statePath, ids[j])
					if(computeLLR) {
						llr <- updateFun$computeLogLikRatio(gr, emit)
						##elementMetadata(gr)$LLR <- updateFun$computeLogLikRatio(gr, emit)
					}
					grl[[j]] <- gr
				}
				break()
			}
		}
		if(returnEmission)  emitArray[, j, ] <- emit
	}
	if(returnEmission) {
		result <- emitArray
	} else {
		grl <- GRangesList(grl)
		names(grl) <- colnames(r)
		result <- unlist(grl)
	}
	##ans <- list(grl=result, paramsB=paramsB, paramsC=paramsC)
	##return(ans)
	return(result)
}

viterbiWrapperG2 <- function(index.samples,
			     cnStates,
			     p.hom=0.05,
			     is.log,
			     limits,
			     normalIndex=3L,
			     nupdates=10,
			     tolerance=5,
			     computeLLR=TRUE,
			     returnEmission=FALSE,
			     verbose=FALSE,
			     grFun,
			     matrixFun,
			     is.snp,
			     anyNP,
			     ...){
	toGRanges <- grFun$toGRanges
	transitionProbs <- grFun$transitionProbs
	tau <- grFun$tau
	b.r <- matrixFun(index.samples)
	g <- b.r[["g"]]
	r <- b.r[["r"]]
	rm(b.r, grFun); gc()
	nc <- ncol(r); nr <- nrow(r); S <- length(cnStates)
	ids <- colnames(r)
	r <- thresholdCopyNumber(r, limits=limits)
	updateFun <- generatorFunG2(r,
				    cnStates=cnStates,
				    normalIndex=normalIndex,
				    tau=tau,
				    limits=limits,
				    p.hom=p.hom,
				    is.log=is.log,
				    computeLLR=computeLLR,
				    verbose=verbose,
				    transitionProbs=transitionProbs)
	statePath <- matrix(NA, nrow(r), nc)
	grl <- vector("list", nc)
	emitG <- gtEmission(object=g,
			    is.snp=is.snp,
			    rohIndex=c(2L, 4L),
			    prGtHom=c(2/3, 0.99, 0.7, 0.99, 1/2, 2/5),
			    S=length(cnStates),
			    log.it=FALSE)
	if(returnEmission) emitArray <- array(NA, dim=c(nrow(r), nc, length(cnStates)))
	## Things that change with iterator i:
	##  --  first two moments
	##  --  mixture probabilities
	if(anyNP) emitB <- matrix(1, nrow(r), length(cnStates))
	for(j in seq_len(nc)){
		paramsC <- updateFun$initialCnParams(j)
		cn <- r[, j]
		sameStatePath <- 0
		llr <- loglik <- rep(0, nupdates)
		for(i in seq_len(nupdates)){
			paramsCprev <- paramsC
			emitC <- updateFun$updateCnEmission(cn, paramsC, j)
			emit <- emitG[,j,]*emitC
			rm(emitC)
			vit.res <- updateFun$fitViterbi(emit)
			if(i > 1) sameStatePath <- if(identical(statePath, vit.res[["statePath"]])) sameStatePath+1 else 0
			statePath <- vit.res[["statePath"]]
			h <- matrix(vit.res[["fv"]]*vit.res[["bv"]],
				    nr, S)
			paramsC <- updateFun$updateCnParams(cn, paramsC, h, j)
			if(verbose) print(paramsC)
			## can we compare the scale factor across runs?
			## Equation 103, Rabiner
			loglik[i] <- sum(log(vit.res[["sf"]]))
			if(verbose) message(loglik[i])
			if(i == 1) next()
			llr[i] <- loglik[i]-loglik[i-1]
			cond <- (abs(llr[i]) < tolerance) || (sameStatePath > 1)
			if(cond | i == nupdates) {
				statePath <- vit.res[["statePath"]]
				gr <- toGRanges(statePath, ids[j])
				if(computeLLR) {
					llr <- updateFun$computeLogLikRatio(gr, emit)
				}
				grl[[j]] <- gr
				break()
			}
		}
		if(returnEmission)  emitArray[, j, ] <- emit
	}
	if(returnEmission) {
		result <- emitArray
	} else {
		grl <- GRangesList(grl)
		names(grl) <- colnames(r)
		result <- grl
	}
	return(unlist(result))
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

viterbiForSnpSet2 <- function(index.samples,
			      grFun,
			      matrixFun,
			      S=2L,
			      rohIndex,
			      normalIndex=1L,
			      computeLLR=TRUE,
			      verbose=FALSE,
			      ...){
	toGRanges <- grFun$toGRanges
	transitionProbs <- grFun$transitionProbs
	tau <- grFun$tau
	g <- matrixFun(index.samples)
	nc <- ncol(g); nr <- nrow(g)
	ids <- colnames(g)
	statePath <- matrix(NA, nrow(g), nc)
	emitG <- gtEmissionFromMatrix(object=g,
				      rohIndex=rohIndex,
				      is.snp=rep(TRUE, nr),
				      S=S,
				      log.it=FALSE,
				      ...)
	updateFun <- generatorFunSnpSet(g=g,
					normalIndex=1L,
					tau=tau,
					computeLLR=computeLLR,
					verbose=verbose,
					transitionProbs=transitionProbs,
					S=S)
	grl <- vector("list", nc)
	J <- ncol(g)
	for(j in seq_len(J)){
		vit.res <- updateFun$fitViterbi(emitG[,j, ])
		statePath <- vit.res[["statePath"]]
		gr <- toGRanges(statePath, j)
		if(computeLLR) elementMetadata(gr)$LLR <- updateFun$computeLogLikRatio(gr, emitG[,j,])
		grl[[j]] <- gr
	}
	grl <- GRangesList(grl)
	names(grl) <- colnames(g)
	unlist(grl)
}

viterbiForSnpSetIce <- function(index.samples,
				grFun,
				matrixFun,
				S=2L,
				normalIndex=1L,
				rohIndex=normalIndex+1L,
				computeLLR=TRUE,
				verbose=FALSE,
				annotationPkg,
				...){
	toGRanges <- grFun$toGRanges
	transitionProbs <- grFun$transitionProbs
	tau <- grFun$tau
	res <- matrixFun(index.samples)
	g <- res[["g"]]
	gt.conf <- res[["gt.conf"]]
	rm(res);gc()
	nc <- ncol(g); nr <- nrow(g)
	ids <- colnames(g)
	statePath <- matrix(NA, nrow(g), nc)
	emitG <- gtEmissionFromMatrix(object=g,
				      gt.conf=gt.conf,
				      ICE=TRUE,
				      rohIndex=rohIndex,
				      is.snp=rep(TRUE, nr),
				      S=S,
				      log.it=FALSE,
				      annotationPkg=annotationPkg,
				      ...)
	updateFun <- generatorFunSnpSet(g=g,
					normalIndex=1L,
					tau=tau,
					computeLLR=computeLLR,
					verbose=verbose,
					transitionProbs=transitionProbs,
					S=S)
	grl <- vector("list", nc)
	J <- ncol(g)
	for(j in seq_len(J)){
		vit.res <- updateFun$fitViterbi(emitG[,j, ])
		statePath <- vit.res[["statePath"]]
		gr <- toGRanges(statePath, j)
		if(computeLLR) elementMetadata(gr)$LLR <- updateFun$computeLogLikRatio(gr, emitG[,j,])
		grl[[j]] <- gr
	}
	grl <- GRangesList(grl)
	names(grl) <- colnames(g)
	unlist(grl)
}
