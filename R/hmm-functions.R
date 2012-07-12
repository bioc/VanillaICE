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
##	if(missing(cnStates)){
##		cnStates <- if(is.log) logCnStates(med.cn) else copyNumberStates(med.cn)
##	}
	limits <- copyNumberLimits(is.log)
##	if(!all.equal(cnStates[normalIndex], med.cn, tolerance=0.2)==TRUE){
##		warning("The median copy number is not near cnStates[3]. The initial values for the mean copy number states are not well specified.")
##	}
	## following is expensive for large data
	if(!missing(sampleIds)) {
		J <- match(sampleIds, sampleNames(object))
	} else J <- seq_len(ncol(object))
	object <- chromosomePositionOrder(object)
	arm <- oligoClasses:::.getArm(chromosome(object), position(object), genomeBuild(object))
	arm <- factor(arm, unique(arm))
	index <- split(seq_len(nrow(object)), arm)
	l <- elementLengths(index)
	if(any(l < 10)) index <- index[l >= 10]
	if(!length(index) > 0) stop("Chromosome arms must have at least 10 markers")
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

hmmBeadStudioSet <- function(object,
			     cnStates=c(-2, -0.4, 0, 0, 0.5, 1),
			     normalIndex=3L,
			     prOutlierCN=0.01,
			     p.hom=0.05,
			     TAUP=1e8,
			     is.log=TRUE,
			     initialProb=rep(1/length(cnStates), length(cnStates)),
			     center=TRUE,
			     nupdates=10,
			     tolerance=5,
			     sampleIds, computeLLR=TRUE, ...){
	## need to check before any subsetting
	if(!missing(cnStates)) {
		if(length(cnStates) != 6) stop("The length of  cnStates should be 6, indicating the initial values for the mean copy number for the states homozygous deletion, hemizygous deletion, normal copy number, normal copy number and no heterozygosity, single-copy amplification, and two-copy amplification.")
		if(cnStates[normalIndex] != cnStates[normalIndex+1])
			stop("the copy number for the normal state and the copy number for regions of homozogyosity should be the same")
	}
	isff <- is(copyNumber(object), "ff")
	##if(missing(is.log)) stop("must specify is.log (TRUE if copy number is on the log-scale)")
	marker.index <- validChromosomeIndex(object)
	gd <- featureData(object)[marker.index, ]
	io <- order(chromosome(gd), position(gd))
	marker.index <- marker.index[io]
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
	##if(missing(is.log)) is.log <- guessCnScale(cnsamp)
	##	if(missing(cnStates)){
	##		cnStates <- if(is.log) logCnStates(med.cn, is.ratio=TRUE) else copyNumberStates(med.cn)
	##	}
	limits <- copyNumberLimits(is.log)
##	if(!all.equal(cnStates[normalIndex], med.cn, tolerance=0.5)==TRUE){
##		message("The median copy number is not near cnStates[normalIndex]. This is reasonable if ")
##	}
	arm <- oligoClasses:::.getArm(chromosome(object)[marker.index], position(object)[marker.index], genomeBuild(object))
	arm <- factor(arm, unique(arm))
	index <- split(marker.index, arm)
	l <- elementLengths(index)
	if(any(l < 10)) index <- index[l >= 10]
	if(!length(index) > 0) stop("Chromosome arms must have at least 10 markers")
	if(!getDoParRegistered()) {
		pkgs <- "VanillaICE"
	}  else {
		pkgs <- c("ff", "VanillaICE")
	}## to avoid warnings
	if(!missing(sampleIds)) {
		J <- match(sampleIds, sampleNames(object))
	} else J <- seq_len(ncol(object))
	i <- NULL
	if(chromosome(object)[1]==23) center <- FALSE  ## do not center for chrom x
	v2fit <- foreach(i=index, .packages="VanillaICE") %do% {
		viterbi2Wrapper(r=copyNumber(object)[i, J, drop=FALSE]/100,
				b=baf(object)[i, J, drop=FALSE]/1000,
				pos=position(object)[i],
				is.snp=isSnp(object)[i],
				chrom=unique(chromosome(object)[i]),
				cnStates=cnStates,
				p.hom=p.hom,
				TAUP=TAUP,
				is.log=is.log,
				center=center,
				limits=limits,
				normalIndex=normalIndex,
				initialProb=initialProb,
				nupdates=nupdates,
				tolerance=tolerance,
				computeLLR=computeLLR, ...)
	}
	## Reorganize so that samples are elements in the list and
	## chromosome arms for the same sample comprise a single
	## GRanges object.
	names(v2fit) <- names(index)
	rd <- stackGRangesList(v2fit, genomeBuild(object))
	return(rd)
}

hmmBeadStudioSetList <- function(object, sampleIds, ...){
	message("Fitting HMM to each chromosome")
	if(isPackageLoaded("ff")) pkgs <- c("ff", "Biobase", "VanillaICE") else pkgs <- c("Biobase", "VanillaICE")
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	##if(ncol(object[[1]]) < ocSamples()){
	if(length(sampleIds) < ocSamples()){
		rdl <- foreach(obj=object, .packages=pkgs) %dopar% {
			hmmBeadStudioSet(object=obj, sampleIds=sampleIds, ...)
		}
		rd <- stackGRangesList(rdl, genomeBuild(object))
		names(rd) <- sampleIds
	} else {
		index <- match(sampleIds, sampleNames(object))
		sample.index <- splitIndicesByLength(index, ocSamples())
		## create a single stream of tasks that can be executed in parallel
		obj <- NULL; j <- NULL
		rdl <- foreach(j=sample.index, .packages=pkgs) %:%
			foreach(obj=object) %dopar%{
				hmmBeadStudioSet(object=obj, sampleIds=sampleNames(object)[j], ...)
			}
		rd <- foreach(ranges=rdl) %do% stackGRangesList(rdl, genomeBuild(object))
##		rdl2 <- foreach(ranges=rdl) %do% stack(RangedDataList(ranges))
##		rdl2 <- foreach(ranges=rdl2) %do% ranges[, -ncol(ranges)]
##		rdl3 <- stack(RangedDataList(rdl2))
##		rd <- rdl3[, -ncol(rdl3)]
	}
	return(rd)
}

hmmSnpSet <- function(object,
		      ICE=FALSE,
		      chromosome=1:22,
		      normalIndex=1L,
		      rohIndex=normalIndex+1L,
		      S=2L,
		      prGtHom=c(.95, 0.05),
		      prGtMis=rep(1/S,S),
		      prHetCalledHom=0.001,
		      prHetCalledHet=0.995,
		      prHomInNormal=0.8,
		      prHomInRoh=0.999, TAUP=1e8, ...){ ## additional arguments to Viterbi construction
	object <- object[chromosome(object) %in% chromosome, ]
	marker.index <- validChromosomeIndex(object)
	object <- object[marker.index, ]
	object <- chromosomePositionOrder(object)
	arm <- oligoClasses:::.getArm(chromosome(object), position(object), genomeBuild(object))
	index <- split(seq_len(nrow(object)), arm)
	l <- elementLengths(index)
	if(any(l < 10)) index <- index[l >= 10]
	if(!length(index) > 0) stop("Chromosome arms must have at least 10 markers")
	chr <- sapply(index, function(i, chr) unique(chr[i]), chr=chromosome(object))
	if(ICE) checkAnnotationForICE(object)
	chrom <- i <- NULL
	rd <- foreach(i=index, chrom=chr) %do% {
		viterbiForSnpSet(gt=calls(object)[i, , drop=FALSE],
				 is.snp=isSnp(object)[i],
				 pos=position(object)[i],
				 gt.conf=confs(object)[i, , drop=FALSE],
				 chrom=chrom,
				 ICE=ICE,
				 annotationPkg=annotation(object),
				 prGtHom=prGtHom,
				 prGtMis=prGtMis,
				 prHetCalledHom=prHetCalledHom,
				 prHetCalledHet=prHetCalledHet,
				 prHomInNormal=prHomInNormal,
				 prHomInRoh=prHomInRoh,
				 TAUP=TAUP,...)
	}
	rd <- stackGRangesList(rd, genomeBuild(object))
	return(rd)
}

hmmOligoSetList <- function(object, sampleIds, ...){
	message("Fitting HMM to each chromosome")
	if(isPackageLoaded("ff")) pkgs <- c("ff", "Biobase", "VanillaICE") else pkgs <- c("Biobase", "VanillaICE")
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	if(ncol(object[[1]]) < ocSamples()){
		rdl <- foreach(obj=object, .packages=pkgs) %dopar% {
			hmmOligoSnpSet(object=obj, sampleIds=sampleIds, ...)
		}
		rd <- stackRangedData(rdl)
	} else {
		sample.index <- splitIndicesByLength(seq_along(sampleIds), ocSamples())
		## create a single stream of tasks that can be executed in parallel
		obj <- NULL; j <- NULL
		rdl <- foreach(j=sample.index, .packages=pkgs) %:%
			foreach(obj=object) %dopar%{
				hmmOligoSnpSet(object=obj, sampleIds=sampleNames(object)[j], ...)
			}
		rdl2 <- foreach(ranges=rdl) %do% stack(RangedDataList(ranges))
		rdl2 <- foreach(ranges=rdl2) %do% ranges[, -ncol(ranges)]
		rdl3 <- stack(RangedDataList(rdl2))
		rd <- rdl3[, -ncol(rdl3)]
	}
	return(rd)
}
