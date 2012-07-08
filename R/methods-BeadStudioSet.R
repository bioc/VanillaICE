setMethod("calls", signature(object="BeadStudioSet"), function(object) assayDataElement(object, "call"))

##setMethod("isInteger", signature(object="BeadStudioSet"), function(object) object@is.integer)

setAs("BeadStudioSet", "data.frame",
      function(from, to){
	      ##isint <- isInteger(from)
	      isint <- TRUE
	      logr.ratio <- as.numeric(lrr(from))
	      if(isint) logr.ratio <- logr.ratio/100
	      bf <- as.numeric(assayDataElement(from, "baf"))
	      if(isint) bf <- bf/1000
	      gt.present <- "call" %in% ls(assayData(from))
	      if(gt.present){
		      gt <- as.integer(calls(from))
	      }
	      x <- rep(position(from)/1e6, ncol(from))
	      ##x <- rep(position(object)[marker.index], 4)/1e6
	      is.snp <- rep(isSnp(from), ncol(from))
	      id <- rep(sampleNames(from), each=nrow(from))
	      if(!gt.present){
		      df <- data.frame(x=x,
				       lrr=logr.ratio,
				       baf=bf,
				       id=id,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE)
	      } else {
		      df <- data.frame(x=x,
				       lrr=logr.ratio,
				       gt=gt,
				       baf=bf,
				       id=id,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE)
	      }
	      return(df)
      })

BeadStudioSet <- function(filenames,
			  lrr.colname="Log.R.Ratio",
			  baf.colname="B.Allele",
			  sep="\t",
			  header=TRUE,
			  colClasses,
			  genome=c("hg19", "hg18"),
			  annotationPkg, chromosome=1:22, ...){
	genome <- match.arg(genome)
	if(missing(colClasses))
		colClasses <- getColClasses(filenames[1], lrr.colname=lrr.colname, baf.colname=baf.colname)
	features <- keyOffFirstFile(filename=filenames[1],
				    cdfname=annotationPkg,
				    genome=genome,
				    colClasses=colClasses,
				    lrr.colname=lrr.colname,
				    baf.colname=baf.colname, ...)
	dat <- read.bsfiles(filenames=filenames,
			    lrr.colname=lrr.colname,
			    baf.colname=baf.colname,
			    row.names=1,
			    sep=sep,
			    drop=FALSE, colClasses=colClasses,
			    nrows=nrow(features)+5000, ...) ## there are about 3-4k markers not in the annotation file
	features <- features[features[, "chrom"] %in% chromosome, ]
	dat <- dat[features[, "index"], , , drop=FALSE]
	r <- as.matrix(dat[,1,])
	b <- as.matrix(dat[,2,])
	colnames(r) <- colnames(b) <- basename(filenames)
	fd <- GenomeAnnotatedDataFrameFrom(r, annotationPkg=annotationPkg, genome=genome)
	new("BeadStudioSet",
	    lrr=integerMatrix(as.matrix(dat[,1, ]), 100),
	    baf=integerMatrix(as.matrix(dat[,2, ]), 1000),
	    annotation=annotationPkg,
	    featureData=fd)
}


hmmBeadStudioSet <- function(object,
			     cnStates,
			     normalIndex=3L,
			     prOutlierCN=0.01,
			     p.hom=0.05,
			     TAUP=1e8,
			     is.log,
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
	if(missing(is.log)) is.log <- guessCnScale(cnsamp)
	if(missing(cnStates)){
		cnStates <- if(is.log) logCnStates(med.cn, is.ratio=TRUE) else copyNumberStates(med.cn)
	}
	limits <- copyNumberLimits(is.log)
	if(!all.equal(cnStates[normalIndex], med.cn, tolerance=0.5)==TRUE){
		warning("The median copy number is not near cnStates[3]. The initial values for the mean copy number states are not well specified.")
	}
	arm <- .getArm(chromosome(object)[marker.index], position(object)[marker.index], genomeBuild(object))
	arm <- factor(arm, unique(arm))
	index <- split(marker.index, arm)
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


##hmmBeadStudioSetCancer <- function(object,
##			     cnStates,
##			     normalIndex=3L,
##			     prOutlierCN=0.01,
##			     p.hom=0.05,
##			     TAUP=1e8,
##			     is.log,
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
##	autosomeIndex <- intersect(marker.index, sample(which(chromosome(object) < 23), ns))
##	cnsamp <- copyNumber(object)[sample(which(chromosome(object) < 23), min(10e3, nrow(object))), sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100
##	med.cn <- median(cnsamp, na.rm=TRUE)
##	if(missing(is.log)) is.log <- guessCnScale(cnsamp)
##	if(missing(cnStates)){
##		cnStates <- if(is.log) logCnStates(med.cn) else copyNumberStates(med.cn)
##	}
##	limits <- copyNumberLimits(is.log)
##	if(!all.equal(cnStates[normalIndex], med.cn, tolerance=0.2)==TRUE){
##		warning("The median copy number is not near cnStates[3]. The initial values for the mean copy number states are not well specified.")
##	}
##	limits <- copyNumberLimits(is.log)
##	arm <- .getArm(chromosome(object)[marker.index], position(object)[marker.index])
##	arm <- factor(arm, unique(arm))
##	index <- split(marker.index, arm)
##	if(!getDoParRegistered()) {
##		pkgs <- "VanillaICE"
##	}  else {
##		pkgs <- c("ff", "VanillaICE")
##	}## to avoid warnings
##	if(!missing(sampleIds)) {
##		J <- match(sampleIds, sampleNames(object))
##	} else J <- seq_len(ncol(object))
##	i <- NULL
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
##				center=FALSE,
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
##	names(v2fit) <- names(index)
##	rd <- stackGRangesList(v2fit, genomeBuild(object))
##	return(rd)
##}



