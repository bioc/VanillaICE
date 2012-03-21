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
			  universe=c("", "hg18", "hg19"),
			  annotationPkg, chromosome=1:22, ...){
	if(universe != "")
		universe <- match.arg(universe)
	if(missing(colClasses))
		colClasses <- getColClasses(filenames[1], lrr.colname=lrr.colname, baf.colname=baf.colname)
	features <- keyOffFirstFile(filename=filenames[1], cdfname=annotationPkg, universe=universe, colClasses=colClasses,
				    lrr.colname=lrr.colname, baf.colname=baf.colname, ...)
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
	fd <- GenomeAnnotatedDataFrameFrom(r, annotationPkg=annotationPkg,
					   universe=universe)
	new("BeadStudioSet",
	    lrr=as.matrix(dat[,1, ]),
	    baf=as.matrix(dat[,2, ]),
	    annotation=annotationPkg,
	    featureData=fd)
}


hmmBeadStudioSet <- function(object,
			     cnStates=logCnStates(),
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
	## need to check before any subsetting
	isff <- is(copyNumber(object), "ff")
	##if(missing(is.log)) stop("must specify is.log (TRUE if copy number is on the log-scale)")
	if(cnStates[normalIndex] != cnStates[rohIndex]){
		stop("the copy number for the normal state and the copy number for regions of homozogyosity should be the same")
	}
	marker.index <- validChromosomeIndex(object)
	gd <- featureData(object)[marker.index, ]
	io <- order(chromosome(gd), position(gd))
	marker.index <- marker.index[io]
	##b <- r <- matrix(NA, length(marker.index), ncol(object))
	## This is expensive.
	##	if(!identical(marker.index, seq_len(nrow(object)))){
	##		object <- object[marker.index, ]
	##	}
	ns <- min(25e3, nrow(object))
	autosomeIndex <- intersect(marker.index, sample(which(chromosome(object) < 23), ns))
	if(missing(is.log)) is.log <- guessCnScale(copyNumber(object)[autosomeIndex, sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE]/100)
	if(!is.log){
		##  for BeadStudioSet assumes log scale.  Need to exponentiate default values
		meds <- median(copyNumber(object)[chromosome(object) < 22, 1]/100, na.rm=TRUE)
		near2 <- all.equal(meds, 2, tolerance=0.2)
		if(!is(near2, "logical")) stop("The median copy number of a random sample of the autosomal markers (after division by 100) is not within 2 +/- 0.2.  Please specify cnStates and is.log for this data")
		if(near2){
			cnStates <- c(0.5, 1.5, 2, 2, 2.5, 3)
		}
	}
	limits <- copyNumberLimits(is.log)
	## could be very expensive and would remove replace
	## the ff object with a matrix...
	##isint <- isInteger(object)
	isint <- TRUE
##	if(isff){
##		## We do not want to write into the ff objects.  Leave intact.
##		## do no more than 100 at a time
##		index <- splitIndicesByLength(seq_len(ncol(object)), 100)
	##		if(isint){
##			for(i in seq_along(index)){
##				j <- index[[i]]
##				r[, j] <- thresholdCopyNumber(copyNumber(object)[marker.index, j]/100, limits)
##				b[, j] <- baf(object)[marker.index, j]/1000
##			}
##		} else {
##			for(i in seq_along(index)){
##				j <- index[[i]]
##				r[, j] <- thresholdCopyNumber(copyNumber(object)[marker.index, j], limits)
##			}
##		}
##	} else {
##		if(isInt){
##			r <- thresholdCopyNumber(copyNumber(object)[marker.index, ], limits)/100
##			b <- baf(object)[marker.index, ]/1000
##		} else {
##			r <- thresholdCopyNumber(copyNumber(object)[marker.index, ], limits)
##			b <- baf(object)[marker.index, ]
##		}
##	}
	## potentially expensive
	## object <- chromosomePositionOrder(object)
	##io <- order(chromosome(object), position(object))
	arm <- .getArm(chromosome(object)[marker.index], position(object)[marker.index])
	##io <- io[marker.index]
	##arm <- arm[marker.index]
	##index <- split(io, arm)
	index <- split(marker.index, arm)
	if(length(index[[1]]) < 1000){
		## better to combine arms for estimates of mean/sd
		index <- list(marker.index)
	}
	if(!getDoParRegistered()) {
		registerDoSEQ()
		pkgs <- "VanillaICE"
	}  else {
		pkgs <- c("ff", "VanillaICE")
	}## to avoid warnings
	## parallelization should be in the calling hmm method, not here.
	v2fit <- foreach(i=index, .packages="VanillaICE") %do% {
		viterbi2Wrapper(r=copyNumber(object)[i, , drop=FALSE]/100,
				b=baf(object)[i, , drop=FALSE]/1000,
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



