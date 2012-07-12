setMethod("calls", signature(object="BeadStudioSet"), function(object) assayDataElement(object, "call"))

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




