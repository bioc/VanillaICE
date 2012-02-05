setMethod("hmm2", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, use.baf=FALSE, k=5, ...){
		  verbose <- hmm.params[["verbose"]] > 0
		  log.beta.cn <- cnEmission(object,
					    cnStates=hmm.params[["copynumberStates"]],
					    k=k,
					    is.log=hmm.params[["is.log"]],
					    is.snp=isSnp(object),
					    normalIndex=hmm.params[["normalIndex"]],...)
		  if(use.baf){
			  if(!"baf" %in% ls(assayData(object))){
				  stop("use.baf is true, but baf not in assayData.  See calculateRBaf.")
			  }
			  log.beta.gt <- bafEmission(object, ...)
		  } else {
			  log.beta.gt <- gtEmission(object, hmm.params, ...)
		  }
		  ##log.emission <- emit(object, hmm.params)
		  log.beta <- log.beta.cn+log.beta.gt
		  ##dimnames(log.beta) <- NULL
		  dimnames(log.beta) <- list(NULL,
					     sampleNames(object),
					     hmm.params$states)
		  res <- viterbi(object, hmm.params, log.E=log.beta)
		  return(res)
})


hmmOligoSnpSet <- function(object,
			   cnStates=c(0, 1, 2, 2, 3, 4),
			   normalIndex=3L,
			   rohIndex=normalIndex+1L,
			   prOutlierCN=0.01,
			   prOutlierBAF=1e-3,
			   p.hom=0.05,
			   TAUP=1e8,
			   is.log,
			   initialProb=rep(1/length(cnStates), length(cnStates)),
			   center=TRUE,
			   reestimation=TRUE,
			   nupdates=10,
			   tolerance=1, ...){
	##if(missing(is.log)) stop("must specify is.log (TRUE if copy number is on the log-scale)")
	if(cnStates[normalIndex] != cnStates[rohIndex]){
		stop("the copy number for the normal state and the copy number for regions of homozogyosity should be the same")
	}
	marker.index <- validChromosomeIndex(object)
	object <- object[marker.index, ]
	if(missing(is.log)) is.log <- guessCnScale(copyNumber(object)[sample(which(chromosome(object) < 23), min(10e3, nrow(object))), sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE])
	if(is.log) cnStates <- logCnStates()
	limits <- copyNumberLimits(is.log)
	copyNumber(object) <- thresholdCopyNumber(copyNumber(object), limits)
	object <- chromosomePositionOrder(object)
	arm <- .getArm(chromosome(object), position(object))
	index <- split(seq_len(nrow(object)), arm)
	v2fit <- foreach(i=index, .packages="VanillaICE", .inorder=TRUE) %do% {
		viterbi2Wrapper(r=copyNumber(object)[i, , drop=FALSE],
				gt=calls(object)[i, , drop=FALSE],
				pos=position(object)[i],
				is.snp=isSnp(object)[i],
				chrom=unique(chromosome(object)[i]),
				cnStates=cnStates,
				prOutlierBAF=prOutlierBAF,
				p.hom=p.hom,
				TAUP=TAUP,
				is.log=is.log,
				center=center,
				reestimation=reestimation,
				limits=limits,
				normalIndex=normalIndex,
				rohIndex=rohIndex,
				initialProb=initialProb,
				nupdates=nupdates,
				tolerance=tolerance)
	}
	rd <- stackRangedData(lapply(v2fit, "[[", 1))
	if(FALSE){
		loglik <- lapply(v2fit, "[[", 2)
		names(loglik) <- paste("chr", names(index))
	}
	return(rd)
}


setMethod("sd", signature(x="oligoSnpSet"),
	  function(x, na.rm=FALSE){
		  getSds(x, na.rm=TRUE)
	   })

setMethod("cnEmission", signature(object="oligoSnpSet"),
	  function(object, stdev, k=5, cnStates, is.log, is.snp, normalIndex, ...){
		  ##fn <- featureNames(object)
		  is.ordered <- checkOrder(object)
		  stopifnot(is.ordered)
		  CN <- copyNumber(object)
		  sds <- sd(object)
		  emit <- cnEmission(object=CN,
				     stdev=sds,
				     k=k,
				     cnStates=cnStates,
				     is.log=is.log,
				     is.snp=is.snp,
				     normalIndex=normalIndex, ...)
		  return(emit)
	  })

setMethod("sd", signature(x="oligoSnpSet"),
	  function(x, na.rm=FALSE){
		  getSds(x, na.rm=TRUE)
	   })



