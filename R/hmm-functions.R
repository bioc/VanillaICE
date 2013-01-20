##hmmOligoSnpSet <- function(object,
##			   cnStates= c(0, 1, 2, 2, 3, 4),
##			   normalIndex=3L,
##			   prOutlierCN=0.01,
##			   TAUP=1e10,
##			   is.log,
##			   initialProb=rep(1/length(cnStates), length(cnStates)),
##			   center=TRUE,
##			   nupdates=10,
##			   tolerance=5,
##			   sampleIds,
##			   computeLLR=TRUE, ...){
##	##if(missing(is.log)) stop("must specify is.log (TRUE if copy number is on the log-scale)")
##	if(!missing(cnStates)){
##		if(length(cnStates) != 6) stop("The length of  cnStates should be 6, indicating the initial values for the mean copy number for the states homozygous deletion, hemizygous deletion, normal copy number, normal copy number and no heterozygosity, single-copy amplification, and two-copy amplification.")
##		if(cnStates[normalIndex] != cnStates[4]) stop("the copy number for the normal state and the copy number for regions of homozogyosity should be the same")
##	}
##	marker.index <- validChromosomeIndex(object)
##	##object <- object[marker.index, ]
##	ns <- min(25e3, nrow(object))
##	if(any(chromosome(object) < 23)){
##		autosomeIndex <- intersect(marker.index, sample(which(chromosome(object) < 23), min(10e3,ns)))
##		cnsamp <- copyNumber(object)[autosomeIndex, sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100
##	} else {
##		if(!any(chromosome(object)==23))
##			return(NULL)
##		cnsamp <- copyNumber(object)[sample(seq_len(nrow(object)), min(10e3, nrow(object))), sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100
##	}
##	med.cn <- median(cnsamp, na.rm=TRUE)
##	if(missing(is.log)) is.log <- guessCnScale(cnsamp)
####	if(missing(cnStates)){
####		cnStates <- if(is.log) logCnStates(med.cn) else copyNumberStates(med.cn)
####	}
##	limits <- copyNumberLimits(is.log)
####	if(!all.equal(cnStates[normalIndex], med.cn, tolerance=0.2)==TRUE){
####		warning("The median copy number is not near cnStates[3]. The initial values for the mean copy number states are not well specified.")
####	}
##	## following is expensive for large data
##	if(!missing(sampleIds)) {
##		J <- match(sampleIds, sampleNames(object))
##	} else J <- seq_len(ncol(object))
##	object <- chromosomePositionOrder(object)
##	arm <- oligoClasses:::.getArm(chromosome(object), position(object), genomeBuild(object))
##	arm <- factor(arm, unique(arm))
##	index <- split(seq_len(nrow(object)), arm)
##	l <- elementLengths(index)
##	if(any(l < 10)) index <- index[l >= 10]
##	if(!length(index) > 0) stop("Chromosome arms must have at least 10 markers")
##	if(chromosome(object)[1]==23) center <- FALSE  ## do not center for chrom x
##	i <- NULL
##	v2fit <- foreach(i=index, .packages="VanillaICE") %do% {
##		viterbiWrapperG(r=copyNumber(object)[i, J, drop=FALSE]/100,
##				gt=calls(object)[i, J, drop=FALSE],
##				pos=position(object)[i],
##				is.snp=isSnp(object)[i],
##				chrom=unique(chromosome(object)[i]),
##				cnStates=cnStates,
##				TAUP=TAUP,
##				is.log=is.log,
##				center=center,
##				limits=limits,
##				normalIndex=normalIndex,
##				initialProb=initialProb,
##				nupdates=nupdates,
##				tolerance=tolerance,
##				computeLLR=computeLLR,...)
##	}
##	names(v2fit) <- names(index)
##	rd <- stackGRangesList(v2fit, genomeBuild(object))
##	return(rd)
##}

## we use genotypes and not the BAFs
hmmOligoSnpSet2 <- function(object, sampleIds, TAUP=1e10, tauMAX,
			    cnStates=c(0, 1, 2, 2, 3, 4),
			    is.log=FALSE,
			    ...){
	pkgs <- c("GenomicRanges", "VanillaICE", "oligoClasses", "matrixStats")
	if(isPackageLoaded("ff")) pkgs <- c("ff", pkgs)
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	index.samples <- match(sampleIds, sampleNames(object))
	chunks <- splitIndicesByLength(index.samples, ocSamples())
	r <- copyNumber(object)
	g <- calls(object)
	pos <- position(object)
	chr <- chromosome(object)
	grFun <- generatorGRanges(chr, pos, genomeBuild(object), TAUP=TAUP, tauMAX=tauMAX)
	is.snp <- isSnp(object)
	snp.index <- which(is.snp)
	anyNP <- any(!is.snp)
	##if(anyNP) np.index <- which(!is.snp) else np.index <- integer()
	center <- TRUE
	pkgs <- c("oligoClasses", "VanillaICE")
	isff <- is(g, "ff")
	if(isff) pkgs <- c("ff", pkgs)
	if(is.log) limits <- c(-4,3) else limits <- c(0, 5)
	matrixFun <- generatorViterbiCnGt(r, g, chr, center=TRUE,
					  snp.index=snp.index,
					  anyNP=anyNP,
					  cnStates=cnStates)
	rm(pos, chr, r, g); gc()
	i <- NULL
	results <- foreach(i=chunks, .packages=pkgs) %dopar% {
		viterbiWrapperG2(index.samples=i,
				 is.snp=is.snp,
				 anyNP=anyNP,
				 is.log=is.log,
				 limits=limits,
				 cnStates=cnStates,
				 grFun=grFun,
				 matrixFun=matrixFun, ...)
	}
	grl <- GRangesList(results)
	metadata(grl) <- list(genome=genomeBuild(object))
	grl
}

##hmmBeadStudioSet <- function(object,
##			     cnStates=c(-2, -0.4, 0, 0, 0.4, 1),
##			     normalIndex=3L,
##			     prOutlierCN=0.01,
##			     p.hom=0.05,
##			     TAUP=1e10,
##			     is.log=TRUE,
##			     initialProb=rep(1/length(cnStates), length(cnStates)),
##			     center=TRUE,
##			     nupdates=10,
##			     tolerance=5,
##			     sampleIds, computeLLR=TRUE, ...){
##	## need to check before any subsetting
##	if(!missing(cnStates)) {
##		if(length(cnStates) != 6) stop("The length of  cnStates should be 6, indicating the initial values for the mean copy number for the states homozygous deletion, hemizygous deletion, normal copy number, normal copy number and no heterozygosity, single-copy amplification, and two-copy amplification.")
##		if(cnStates[normalIndex] != cnStates[normalIndex+1])
##			stop("the copy number for the normal state and the copy number for regions of homozogyosity should be the same")
##	}
##	isff <- is(copyNumber(object), "ff")
##	##if(missing(is.log)) stop("must specify is.log (TRUE if copy number is on the log-scale)")
##	marker.index <- validChromosomeIndex(object)
##	gd <- featureData(object)[marker.index, ]
##	io <- order(chromosome(gd), position(gd))
##	marker.index <- marker.index[io]
##	ns <- min(25e3, nrow(object))
##	if(any(chromosome(object) < 23)){
##		autosomeIndex <- intersect(marker.index, sample(which(chromosome(object) < 23), min(10e3,ns)))
##		cnsamp <- copyNumber(object)[autosomeIndex, sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100
##	} else {
##		if(!any(chromosome(object)==23))
##			return(NULL)
##		cnsamp <- copyNumber(object)[sample(seq_len(nrow(object)), min(10e3, nrow(object))), sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100
##	}
##	med.cn <- median(cnsamp, na.rm=TRUE)
##	##if(missing(is.log)) is.log <- guessCnScale(cnsamp)
##	##	if(missing(cnStates)){
##	##		cnStates <- if(is.log) logCnStates(med.cn, is.ratio=TRUE) else copyNumberStates(med.cn)
##	##	}
##	limits <- copyNumberLimits(is.log)
####	if(!all.equal(cnStates[normalIndex], med.cn, tolerance=0.5)==TRUE){
####		message("The median copy number is not near cnStates[normalIndex]. This is reasonable if ")
####	}
##	arm <- oligoClasses:::.getArm(chromosome(object)[marker.index], position(object)[marker.index], genomeBuild(object))
##	arm <- factor(arm, unique(arm))
##	if(any(table(arm) < 1000)){
##		drop.level <- names(table(arm))[table(arm) < 1000]
##		if(length(grep("p", drop.level)) > 0){
##			collapselevel <- gsub("p", "q", drop.level[grep("p", drop.level)])
##			for(i in seq_along(drop.level)){
##				arm[arm == drop.level[i]] <- collapselevel[i]
##			}
##		}
##		if(length(grep("q", drop.level)) > 0){
##			collapselevel <- gsub("q", "p", drop.level[grep("q", drop.level)])
##			for(i in seq_along(drop.level)){
##				arm[arm == drop.level[i]] <- collapselevel[i]
##			}
##		}
##		arm <- factor(arm, unique(arm))
##	}
##	index <- split(marker.index, arm)
##	l <- elementLengths(index)
##	if(any(l < 10)) index <- index[l >= 10]
##	if(!length(index) > 0) stop("Chromosome arms must have at least 10 markers")
##	if(!getDoParRegistered()) {
##		pkgs <- "VanillaICE"
##	}  else {
##		pkgs <- c("ff", "VanillaICE")
##	}## to avoid warnings
##	if(!missing(sampleIds)) {
##		J <- match(sampleIds, sampleNames(object))
##	} else J <- seq_len(ncol(object))
##	i <- NULL
##	if(chromosome(object)[1]==23) center <- FALSE  ## do not center for chrom x
##	## if object contains data for one chromosome, it is inefficient to unnecessarily subset it
##	v2fit <- foreach(i=index, .packages="VanillaICE") %do% {
##		viterbi2Wrapper(r=copyNumber(object)[i, J, drop=FALSE]/100,
##				b=baf(object)[i, J, drop=FALSE]/1000,
##				pos=position(object)[i],
##				is.snp=isSnp(object)[i],
##				chrom=unique(chromosome(object)[i]),
##				cnStates=cnStates,
##				p.hom=p.hom,
##				TAUP=TAUP,
##				is.log=is.log,
##				center=center,
##				limits=limits,
##				normalIndex=normalIndex,
##				initialProb=initialProb,
##				nupdates=nupdates,
##				tolerance=tolerance,
##				computeLLR=computeLLR, ...)
##	}
##	## Reorganize so that samples are elements in the list and
##	## chromosome arms for the same sample comprise a single
##	## GRanges object.
##	grl <- lapply(v2fit, "[[", 1)
##	rd <- stackGRangesList(grl, genomeBuild(object))
##	paramsB <- lapply(v2fit, "[[", 2)
##	paramsC <- lapply(v2fit, "[[", 3)
##	names(paramsC) <- names(paramsB) <- names(index)
##	result <- list(grl=rd, paramsB=paramsB, paramsC=paramsC)
##	##names(v2fit) <- names(index)
##	##rd <- stackGRangesList(v2fit[["grl"]], genomeBuild(object))
##	##return(rd)
##	return(result)
##}



hmmBafLrrSet2 <- function(object, sampleIds, TAUP=1e10, tauMAX,
			  cnStates=c(-2, -0.4, 0, 0, 0.4, 1),
			  is.log=TRUE,
			  ...){
	pkgs <- c("GenomicRanges", "VanillaICE", "oligoClasses", "matrixStats")
	if(isPackageLoaded("ff")) pkgs <- c("ff", pkgs)
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	index.samples <- match(sampleIds, sampleNames(object))
	chunks <- splitIndicesByLength(index.samples, ocSamples())
	r <- copyNumber(object)
	b <- baf(object)
	pos <- position(object)
	chr <- chromosome(object)
	grFun <- generatorGRanges(chr, pos, genomeBuild(object), TAUP=TAUP, tauMAX=tauMAX)
	is.snp <- isSnp(object)
	snp.index <- which(is.snp)
	anyNP <- any(!is.snp)
	##if(anyNP) np.index <- which(!is.snp) else np.index <- integer()
	center <- TRUE
	pkgs <- c("oligoClasses", "VanillaICE")
	isff <- is(r, "ff")
	if(isff) pkgs <- c("ff", pkgs)
	matrixFun <- generatorViterbiBeadStudioSet(r, b, chr, center=TRUE,
						   snp.index=snp.index, anyNP=anyNP,
						   cnStates=cnStates)
	if(is.log) limits <- c(-4,3) else limits <- c(0, 5)
	rm(pos, chr, b, r); gc()
	i <- NULL
	results <- foreach(i=chunks, .packages=pkgs) %dopar% {
		viterbi2Wrapper(index.samples=i,
				snp.index=snp.index,
				anyNP=anyNP,
				is.log=is.log,
				limits=limits,
				cnStates=cnStates,
				grFun=grFun,
				matrixFun=matrixFun, ...)
	}
	grl <- GRangesList(results)
	metadata(grl) <- list(genome=genomeBuild(object))
	grl
}

##fithmm <- function(object, j, sample.size, ...){
##	grl <- list()
##	m <- 1
##	for(k in seq_along(object)){
##		obj <- object[[k]] ## one chromosome
##		## hmmBeadStudioSet has overhead -- inefficient if there are many samples
##		sb <- splitIndicesByLength(j, sample.size)
##		for(i in seq_along(sb)){
##			jj <- sb[[i]]
##			grl[[m]] <- hmmBeadStudioSet(obj[, jj], ...)
##			elementMetadata(grl[[m]])$sample <- sampleNames(obj)[jj]
##			m <- m+1
##		}
##	}
##	gr <- stackGRangesList(grl, genomeBuild(object))
##	return(gr)
##}
##
##hmmBeadStudioSetList <- function(object, sampleIds, ...){
##	message("Fitting HMM to each chromosome")
##	pkgs <- c("GenomicRanges", "Biobase", "VanillaICE", "oligoClasses")
##	if(isPackageLoaded("ff")) pkgs <- c("ff", pkgs)
##	if(missing(sampleIds)) sampleIds <- sampleNames(object)
##	sample.index <- match(sampleIds, sampleNames(object))
##	if(any(is.na(sample.index))) stop(paste("the following sampleIds are not in sampleNames(object): ", sampleIds[is.na(sampleIds)]))
##	chromBatches <- splitIndicesByNode(seq_along(object))
##	sample.size <- ocSamples()
##	grl <- foreach(i=chromBatches, .packages=pkgs) %dopar% VanillaICE:::fithmm(object=object[i], sample.index, sample.size=sample.size, ...)
##	gr <- stackGRangesList(grl, genomeBuild(object))
##	return(gr)
##}
##
##
##
##
##hmmSnpSet <- function(object,
##		      ICE=FALSE,
##		      chromosome=1:22,
##		      normalIndex=1L,
##		      rohIndex=normalIndex+1L,
##		      S=2L,
##		      prGtHom=c(.95, 0.05),
##		      prGtMis=rep(1/S,S),
##		      prHetCalledHom=0.001,
##		      prHetCalledHet=0.995,
##		      prHomInNormal=0.8,
##		      prHomInRoh=0.999, TAUP=1e10, ...){ ## additional arguments to Viterbi construction
##	object <- object[chromosome(object) %in% chromosome, ]
##	marker.index <- validChromosomeIndex(object)
##	object <- object[marker.index, ]
##	object <- chromosomePositionOrder(object)
##	arm <- oligoClasses:::.getArm(chromosome(object), position(object), genomeBuild(object))
##	index <- split(seq_len(nrow(object)), arm)
##	l <- elementLengths(index)
##	if(any(l < 10)) index <- index[l >= 10]
##	if(!length(index) > 0) stop("Chromosome arms must have at least 10 markers")
##	chr <- sapply(index, function(i, chr) unique(chr[i]), chr=chromosome(object))
##	if(ICE) checkAnnotationForICE(object)
##	chrom <- i <- NULL
##	rd <- foreach(i=index, chrom=chr) %do% {
##		viterbiForSnpSet(gt=calls(object)[i, , drop=FALSE],
##				 is.snp=isSnp(object)[i],
##				 pos=position(object)[i],
##				 gt.conf=confs(object)[i, , drop=FALSE],
##				 chrom=chrom,
##				 ICE=ICE,
##				 annotationPkg=annotation(object),
##				 prGtHom=prGtHom,
##				 prGtMis=prGtMis,
##				 prHetCalledHom=prHetCalledHom,
##				 prHetCalledHet=prHetCalledHet,
##				 prHomInNormal=prHomInNormal,
##				 prHomInRoh=prHomInRoh,
##				 TAUP=TAUP,...)
##	}
##	rd <- stackGRangesList(rd, genomeBuild(object))
##	return(rd)
##}

hmmSnpSet2 <- function(object, sampleIds, TAUP=1e10, tauMAX,
		       normalIndex=1L,
		       rohIndex=normalIndex+1L,
		       S=2L,
		       ...){
	pkgs <- c("GenomicRanges", "VanillaICE", "oligoClasses", "matrixStats")
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	index.samples <- match(sampleIds, sampleNames(object))
	chunks <- splitIndicesByLength(index.samples, ocSamples())
	g <- calls(object)
	pos <- position(object)
	chr <- chromosome(object)
	grFun <- generatorGRanges(chr, pos, genomeBuild(object), TAUP=TAUP, tauMAX=tauMAX)
	isff <- is(g, "ff")
	if(isff) pkgs <- c("ff", pkgs)
	matrixFun <- generatorViterbiSnpSet(g, chr, isff)
	rm(pos, chr); gc()
	i <- NULL
	results <- foreach(i=chunks, .packages=pkgs) %dopar% {
		viterbiForSnpSet2(index.samples=i,
				  grFun=grFun,
				  matrixFun=matrixFun,
				  S=S,
				  ...)
	}
	grl <- GRangesList(results)
	metadata(grl) <- list(genome=genomeBuild(object))
	grl
}

hmmSnpSetIce <- function(object, sampleIds, TAUP=1e10, tauMAX,
			 ...){
	pkgs <- c("GenomicRanges", "VanillaICE", "oligoClasses", "matrixStats")
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	index.samples <- match(sampleIds, sampleNames(object))
	chunks <- splitIndicesByLength(index.samples, ocSamples())
	g <- calls(object)
	gt.conf <- snpCallProbability(object)
	pos <- position(object)
	chr <- chromosome(object)
	grFun <- generatorGRanges(chr, pos, genomeBuild(object),
				  TAUP=TAUP, tauMAX=tauMAX)
	isff <- is(g, "ff")
	if(isff) pkgs <- c("ff", pkgs)
	matrixFun <- generatorViterbiSnpSetIce(g, gt.conf, chr, isff)
	rm(pos, chr, g, gt.conf); gc()
	i <- NULL
	results <- foreach(i=chunks, .packages=pkgs) %dopar% {
		viterbiForSnpSetIce(index.samples=i,
				    grFun=grFun,
				    matrixFun=matrixFun,
				    annotationPkg=annotation(object),
				    ...)
	}
	GRangesList(results)
}

##hmmOligoSetList <- function(object, sampleIds, ...){
##	message("Fitting HMM to each chromosome")
##	if(isPackageLoaded("ff")) pkgs <- c("ff", "Biobase", "VanillaICE") else pkgs <- c("Biobase", "VanillaICE")
##	if(missing(sampleIds)) sampleIds <- sampleNames(object)
##	##if(ncol(object[[1]]) < ocSamples()){
##	## hmmBeadStudioSet parallelizes over chromosome arms, but not samples
##	## So, when passing a single chromosome, it will at most do two jobs in parallel.
##	sample.index <- match(sampleIds, sampleNames(object))
##	if(any(is.na(sample.index))) stop(paste("the following sampleIds are not in sampleNames(object): ", sampleIds[is.na(sampleIds)]))
##	chromBatches <- splitIndicesByNode(seq_along(object))
##	is.baf <- "baf" %in% ls(assayData(object))
##	fitoligo <- function(object, j, ...){
##		grl <- list()
##		m <- 1
##		for(k in seq_along(object)){
##			obj <- object[[k]] ## one chromosome
##			## hmmBeadStudioSet has overhead -- inefficient if there are many samples
##			sb <- splitIndicesByLength(j, ocSamples())
##			for(i in seq_along(sb)){
##				jj <- sb[[i]]
##				if(is.baf){
##					grl[[m]] <- hmmBeadStudioSet(obj[, jj], ...)
##				} else grl[[m]] <- hmmOligoSnpSet(obj[, jj], ...)
##				m <- m+1
##			}
##		}
##		gr <- stackGRangesList(grl, genomeBuild(object))
##		return(gr)
##	}
##	rdl <- foreach(i=chromBatches, .packages="VanillaICE") %dopar% fitoligo(object[i], sample.index)
##	rdl <- stackGRangesList(rd, genomeBuild(object))
##	return(rd)
##}





## all chromosomes at once.  parallelize over samples
hmmBafLrrSetList2 <- function(object, sampleIds, TAUP=1e10, tauMAX,
			      cnStates=c(-2, -0.4, 0, 0, 0.4, 1), is.log=TRUE,
			      ...){
	pkgs <- c("GenomicRanges", "VanillaICE", "oligoClasses", "matrixStats")
	if(isPackageLoaded("ff")) pkgs <- c("ff", pkgs)
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	index.samples <- match(sampleIds, sampleNames(object))
	chunks <- splitIndicesByLength(index.samples, ocSamples())
	rlist <- lrr(object)
	blist <- baf(object)
	pos <- unlist(position(object))
	chr <- rep(chromosome(object), elementLengths(object))
	grFun <- generatorGRanges(chr, pos, genomeBuild(object), TAUP=TAUP, tauMAX=tauMAX)
	is.snp <- unlist(lapply(featureDataList(object), isSnp))
	snp.index <- which(is.snp)
	anyNP <- any(!is.snp)
	##if(anyNP) np.index <- which(!is.snp) else np.index <- integer()
	center <- TRUE
	pkgs <- c("oligoClasses", "VanillaICE")
	isff <- is(rlist[[1]], "ff")
	if(isff) pkgs <- c("ff", pkgs)
	matrixFun <- generatorViterbi(rlist, blist, chr, center=TRUE,
				      snp.index=snp.index, anyNP=anyNP)
	rm(pos, chr, blist, rlist); gc()
	i <- NULL
	results <- foreach(i=chunks, .packages=pkgs) %dopar% {
		viterbi2Wrapper(index.samples=i,
				snp.index=snp.index,
				anyNP=anyNP,
				is.log=TRUE,
				limits=c(-4, 3),
				cnStates=cnStates,
				grFun=grFun,
				matrixFun=matrixFun, ...)
	}
	grl <- GRangesList(results)
	metadata(grl) <- list(genome=genomeBuild(object))
	grl
}
