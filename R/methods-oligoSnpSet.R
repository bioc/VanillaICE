hmmOligoSnpSet <- function(object,
			   cnStates= c(0, 1, 2, 2, 3, 4),
			   normalIndex=3L,
			   prOutlierCN=0.01,
			   TAUP=1e8,
			   is.log,
			   initialProb=rep(1/length(cnStates), length(cnStates)),
			   center=TRUE,
			   nupdates=10,
			   tolerance=5,
			   sampleIds,
			   computeLLR=TRUE, ...){
	##if(missing(is.log)) stop("must specify is.log (TRUE if copy number is on the log-scale)")
	if(!missing(cnStates)){
		if(length(cnStates) != 6) stop("The length of  cnStates should be 6, indicating the initial values for the mean copy number for the states homozygous deletion, hemizygous deletion, normal copy number, normal copy number and no heterozygosity, single-copy amplification, and two-copy amplification.")
		if(cnStates[normalIndex] != cnStates[4]) stop("the copy number for the normal state and the copy number for regions of homozogyosity should be the same")
	}
	marker.index <- validChromosomeIndex(object)
	##object <- object[marker.index, ]
	ns <- min(25e3, nrow(object))
	if(any(chromosome(object) < 23)){
		autosomeIndex <- intersect(marker.index, sample(which(chromosome(object) < 23), min(10e3,ns)))
		cnsamp <- copyNumber(object)[autosomeIndex, sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100
	} else {
		if(!any(chromosome(object)==23))
			return(NULL)
		cnsamp <- copyNumber(object)[sample(seq_len(nrow(object)), min(10e3, nrow(object))), sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100
	}
	med.cn <- median(cnsamp, na.rm=TRUE)
	if(missing(is.log)) is.log <- guessCnScale(cnsamp)
	if(missing(cnStates)){
		cnStates <- if(is.log) logCnStates(med.cn) else copyNumberStates(med.cn)
	}
	limits <- copyNumberLimits(is.log)
	if(!all.equal(cnStates[normalIndex], med.cn, tolerance=0.2)==TRUE){
		warning("The median copy number is not near cnStates[3]. The initial values for the mean copy number states are not well specified.")
	}
	## following is expensive for large data
	if(!missing(sampleIds)) {
		J <- match(sampleIds, sampleNames(object))
	} else J <- seq_len(ncol(object))
	object <- chromosomePositionOrder(object)
	arm <- .getArm(chromosome(object), position(object), genomeBuild(object))
	arm <- factor(arm, unique(arm))
	index <- split(seq_len(nrow(object)), arm)
	if(chromosome(object)[1]==23) center <- FALSE  ## do not center for chrom x
	i <- NULL
	v2fit <- foreach(i=index, .packages="VanillaICE", .inorder=TRUE) %do% {
		viterbiWrapperG(r=copyNumber(object)[i, J, drop=FALSE]/100,
				gt=calls(object)[i, J, drop=FALSE],
				pos=position(object)[i],
				is.snp=isSnp(object)[i],
				chrom=unique(chromosome(object)[i]),
				cnStates=cnStates,
				TAUP=TAUP,
				is.log=is.log,
				center=center,
				limits=limits,
				normalIndex=normalIndex,
				initialProb=initialProb,
				nupdates=nupdates,
				tolerance=tolerance,
				computeLLR=computeLLR,...)
	}
	names(v2fit) <- names(index)
	rd <- stackGRangesList(v2fit, genomeBuild(object))
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

##setMethod("baf", signature(object="oligoSnpSet"), function(object) assayDataElement(object, "baf"))
