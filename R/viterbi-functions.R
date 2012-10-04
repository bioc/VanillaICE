

viterbi2Wrapper <- function(r, b,
			    pos, is.snp, cnStates,
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
			    verbose=FALSE,
			    ...){
	nc <- ncol(r); nr <- nrow(r); S <- length(cnStates)
	anyNP <- any(!is.snp)
	if(anyNP) {
		snp.index <- which(is.snp)
		if(length(snp.index) == 0) stop("is.snp is all FALSE")
	}
	r <- thresholdCopyNumber(r, limits=limits)
	updateFun <- generatorFun(r, b, is.snp=is.snp, cnStates=cnStates,
				  normalIndex=normalIndex,  TAUP=TAUP,
				  limits=limits, center=center,
				  prOutlierBAF=prOutlierBAF, p.hom=p.hom,
				  position=pos, is.log=is.log,
				  computeLLR=computeLLR, chrom=chrom, verbose=verbose)
	statePath <- matrix(NA, nrow(r), nc)
	grl <- vector("list", nc)
	if(returnEmission) emitArray <- array(NA, dim=c(nrow(r), nc, length(cnStates)))
	## Things that change with iterator i:
	##  --  first two moments
	##  --  mixture probabilities
	## pseudocode
	if(anyNP) emitB <- matrix(1, nrow(r), length(cnStates))
	for(j in seq_len(nc)){
		paramsB <- updateFun$initialBafParams()
		paramsC <- updateFun$initialCnParams(j)
		## if we subset r and bf here we wouldn't need to do "[" for each update i
		bf <- if(anyNP) b[snp.index, j] else b[, j]
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
			if(i > 1){
				emit <- emitB*emitC
			} else emit <- emitB
			vit.res <- updateFun$fitViterbi(emit)
			if(i > 1) sameStatePath <- if(identical(statePath, vit.res[["statePath"]])) sameStatePath+1 else 0
			statePath <- vit.res[["statePath"]]
			fv <- vit.res[["fv"]]
			bv <- vit.res[["bv"]]
			h <- matrix(fv*bv, nr, S)
			if(anyNP){
				paramsB <- updateFun$updateBafParams(bf, paramsB, h[is.snp, ], j)
			} else paramsB <- updateFun$updateBafParams(bf, paramsB, h, j)
			if(verbose) print(paramsB)
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
				gr <- updateFun$toGRanges(statePath, j)
				if(computeLLR) elementMetadata(gr)$LLR <- updateFun$computeLogLikRatio(gr, emit)
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
	return(result)
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
