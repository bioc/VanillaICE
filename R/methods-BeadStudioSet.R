setMethod("calls", signature(object="BeadStudioSet"), function(object) assayDataElement(object, "call"))

setAs("BeadStudioSet", "data.frame",
      function(from, to){
	      logr.ratio <- as.numeric(lrr(from))
	      bf <- as.numeric(assayDataElement(from, "baf"))
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
			     reestimation=TRUE, nupdates=10,
			     tolerance=1, ...){
	##if(missing(is.log)) stop("must specify is.log (TRUE if copy number is on the log-scale)")
	if(cnStates[normalIndex] != cnStates[rohIndex]){
		stop("the copy number for the normal state and the copy number for regions of homozogyosity should be the same")
	}
	marker.index <- validChromosomeIndex(object)
	object <- object[marker.index, ]
	if(missing(is.log)) is.log <- guessCnScale(copyNumber(object)[sample(which(chromosome(object) < 23), min(10e3, nrow(object))), sample(seq_len(ncol(object)), min(2, ncol(object))), drop=FALSE])
	if(!is.log){
		##  for BeadStudioSet assumes log scale.  Need to exponentiate default values
		meds <- median(copyNumber(object)[chromosome(object) < 22, 1], na.rm=TRUE)
		near2 <- all.equal(meds, 2, tolerance=0.2)
		if(near2){
			cnStates <- c(0.5, 1.5, 2, 2, 2.5, 3)
		} else stop("Please specify cnStates and is.log for this data")
	}
	limits <- copyNumberLimits(is.log)
	copyNumber(object) <- thresholdCopyNumber(copyNumber(object), limits)
	object <- chromosomePositionOrder(object)
	arm <- .getArm(chromosome(object), position(object))
	index <- split(seq_len(nrow(object)), arm)
	##indexChrom <- as.integer(sapply(names(index), function(x) strsplit(x, "[pq]")[[1]][[2]]))
	##index <- index[order(indexChrom)]
	v2fit <- foreach(i=index, .packages="VanillaICE") %do% {
		viterbi2Wrapper(r=copyNumber(object)[i, , drop=FALSE],
				b=baf(object)[i, , drop=FALSE],
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



