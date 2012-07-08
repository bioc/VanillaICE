\name{hmm-methods}
\docType{methods}
\alias{hmm}
\alias{hmm-methods}
\alias{hmm,BeadStudioSet-method}
\alias{hmm,CNSet-method}
\alias{hmm,oligoSnpSet-method}
\alias{hmm,SnpSet2-method}
\alias{hmm,BeadStudioSetList-method}
\alias{hmm,oligoSetList-method}
\title{Hidden Markov Model methods}
\description{
 Hidden Markov Model methods in package \pkg{VanillaICE}
}
\section{Methods}{


  \describe{

    The \code{hmm} method is defined for several classes of containers
    of preprocessed and normalized SNP array data.  The most common
    containers for use with genotyping platforms are the
    \code{BeadStudioSet} and \code{oligoSnpSet} classes.  The primary
    difference between these two containers are the requirements for the
    assay data elements.  A \code{BeadStudioSet} object must have assay
    data elements "lrr" (log R ratios) and "baf" (B allele frequencies).
    As of version 1.18.0, all matrices stored in assay data are assumed
    to be integers.  For copy number and relative copy number, the
    estimates should be scaled by 100.  For B allele frequencies, the
    estimates should be scaled by 1000.  The helper function
    \code{integerMatrix} in the oligoClasses package can be useful for
    the conversion. Genotype calls are optional for the
    \code{BeadStudioSet} object.  As the name implies, the
    \code{BeadStudioSet} container would typically be generated as part
    of a pipeline to process data from Illumina array platforms.  By
    contrast, the \code{oligoSnpSet} object has required assay data
    elements "call" (genotype calls), "callProbability" (genotype
    confidence scores), "copyNumber", and "cnConfidence".  As B allele
    frequencies are perhaps more informative than the genotype calls for
    distinguishing copy number states (particularly amplifications), an
    assay data element named "baf" can be included in the assay data for
    an \code{oligoSnpSet} object.  The presence of a "baf" element in
    the assay data of an \code{oligoSnpSet} has implications on the
    particular HMM fit to identify the CNV boundaries (as discussed
    below).

    \item{\code{signature(object = "BeadStudioSet", ...)}}{

      A hidden Markov model for the \code{BeadStudioSet} class.  The
  assay data are log R ratios and B allele frequencies. See
  \code{hmmBeadStudioSet} for additional arguments that can be passed
  through the \code{...} operator.

}

\item{\code{signature(object = "SnpSet2", ...)}}{ A hidden Markov model
  for the \code{SnpSet} class.  The assay data are diallelic genotype
  calls represented as integers (1=AA, 2=AB, 3=BB).  See
  \code{hmmSnpSet} for additional arguments that can be passed
  through the \code{...} operator.
}

\item{\code{signature(object = "CNSet", ...)}}{

  A hidden Markov model for the \code{CNSet} class.  The \code{CNSet}
  instance is first coerced to an object of class \code{oligoSnpSet}
  containing estimates of total copy number and B allele frequencies.
  See \code{hmmBeadStudioSet} for additional arguments that can be
  passed through the \code{...} operator.  For large data sets, the
  initial coercion to the \code{oligoSnpSet} class can be very expensive
  in terms of I/O and require a large amount of RAM.  Users with large
  data sets may prefer to coerce selected samples (e.g., the set of
  samples belonging to a given batch) to an \code{oligoSnpSet} object,
  and then fit the hmm on the \code{oligoSnpSet} object directly.  This
  approach is illustrated in the \code{crlmmDownstream} vignette.

}

\item{\code{signature(object = "CopyNumberSet", ...)}}{

  A hidden Markov model for the \code{CopyNumberSet} class.  The assay
  data are estimates of total copy number.  This method should not be
  used for arrays with genotype information as the genotypes / B allele
  frequencies are informative for copy number inference.

}

\item{\code{signature(object = "oligoSnpSet", ...)}}{

A hidden Markov model for the \code{oligoSnpSet} class.  If "baf" is
included among the assay data elements, the \code{hmmBeadStudioSet} HMM
is implemented.  Otherwise, the \code{hmmOligoSnpSet} is implemented.

}

\item{\code{signature(object = "oligoSetList", ...)}}{

  The oligoSetList class is a container for genotypes, B allele
  frequencies (optional), and copy number organized by chromosome.
  Each element in the list class contains low-level summaries and
  phenotypic information for a single chromosome. The organization by
  chromosome facilitates parallelization of methods to identify copy
  number alterations.  If B allele frequencies are included, the hmm
  fit to instance of this object is the same as the hmm fit to
  instances of a \code{BeadStudioSetList} object (the function
  \code{hmmBeadStudioSet} is fit to each element in the
  \code{oligoSetList} object).

}

\item{\code{signature(object = "BeadStudioSetList", ...)}}{

	The only difference with \code{oligoSetList} is that the
	assayData for \code{BeadStudioSetList} objects must include B
	allele frequencies (B allele frequencies are optional in the
	\code{oligoSetList} class).  The function
	\code{hmmBeadStudioSet} is fit to each element in the
	\code{BeadStudioSetList} object.

}

}
}

\seealso{
	\code{\linkS4class{oligoSetList}}, \code{\linkS4class{BeadStudioSetList}},
	\code{\link{hmmBeadStudioSet}}, \code{\link{hmmOligoSnpSet}},
  \code{\link{hmmSnpSet}}.  For plotting copy number and B allele
  frequencies, see \code{\link[SNPchip]{xyplotLrrBaf}},
  \code{\link[SNPchip]{xypanelBaf}}.

 }

\examples{
	library(oligoClasses)
	library(IRanges)
	data(oligoSetExample, package="oligoClasses")
	oligoSet <- oligoSet[chromosome(oligoSet) == 1, ]
	hmmResults <- hmm(oligoSet)
	state(hmmResults[[1]])
	##
	## Plotting ranges:
	##
	if(require(SNPchip) && require(IRanges)){
		## Plot the data for the second range with a blue
		## border, and frame the region by 10 Mb on each side
		## of the state boundary.
		##
		res <- hmmResults[[1]]
                elementMetadata(res)$sampleId <- names(hmmResults)
		xyplot(cn~x, oligoSet, range=res[2, ], frame=10e6,
		       panel=xypanel, pch=21, cex=0.3,
		       col.hom="royalblue", fill.hom="royalblue",
		       col.het="red", fill.het="red", xlab="Mb",
		       ylab=expression(log[2]("copy number")))
		##  (Note that the formula cn~x is required at this time)
		##
		## Or, plot each range in its own panel with a frame
		##  of 2e6 bases.  (Again, the formula is a standard format
		##  with cn, x, range, and id the only allowed terms) Because
		##  these are all the ranges from one individual's chromosome,
		##  the ranges are overlapping The range 'in focus' is
		##  demarcated by vertical blue lines
		xyplot(cn~x | range, oligoSet, range=res, frame=2e6,
		       panel=xypanel,
		       pch=21,
		       cex=0.3,
		       scales=list(x="free"),
		       border="blue",
		       col.hom="royalblue",
		       col.het="salmon",
		       col.np="grey",
		       par.strip.text=list(cex=0.6),
		       xlab="Mb",
		       ylab=expression(log[2]("copy number")))
	}
	##---------------------------------------------------------------------------
	## For an oligoSnpSet with B allele frequencies:
	##---------------------------------------------------------------------------
	path <- system.file("extdata", package="VanillaICE")
	load(file.path(path, "oligosetForUnitTest.rda"))
	## copy number estimates in this object are not on the log
	## scale, so specify is.log=FALSE and provide the means for
	## the latent copy number states.  IN addition we also specify
	## an initial value and constraints for the probability that
	## the BAF is an outlier
	fit <- hmm(oligoset, is.log=FALSE, cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2),
		   prOutlierBAF=list(initial=1e-4, max=1e-3, maxROH=1e-5))
	##
	##  For log R ratios, one could simply do
	##  hmm(oligoset, prOutlierBAF=list(initial=1e-4, max=1e-3, maxROH=1e-5))
	##
	if(require(SNPchip)){
		## plotting this data
                ## For plotting copy number and log R ratios for multiple genomic intervals, see xyplotLrrBaf
		fit <- fit[[1]]
		library(IRanges)
		library(Biobase)
		rect2 <- function(object){
			col <- c("red", "red", "white", "grey70", "royalblue", "blue")
			object <- object[state(object) !=3 , ]
			object <- object[order(width(object), decreasing=TRUE), ]
			rect(xleft=start(object)/1e6,xright=end(object)/1e6,
			     ybottom=rep(0.7,length(object)),
			     ytop=rep(1,length(object)),
			     col=col[state(object)],
			     border=col[state(object)])
		}
		par(las=1)
		plot(position(oligoset)/1e6, copyNumber(oligoset)/100,
		     pch=".", col="black",
		     ylim=c(-1, 3), ylab="copy number", xlab="position (Mb)")
		rescale <- function(x, l, u){
			b <- 1/(u-l)
			a <- l*b
			(x+a)/b
		}
		b <- rescale(baf(oligoset)/1000, -1, 0)
		rect2(fit)
		franges <- makeFeatureGRanges(oligoset)
		o <- subjectHits(findOverlaps(fit[4, ], franges))
		points(position(oligoset)/1e6, b, pch=".", col="royalblue")
		axis(side=4, at =c(-1, -0.5, 0), labels=c(0, 0.5, 1), col="blue")
		text(10, 0.1, "BAF", col="blue")
	}

	##---------------------------------------------------------------------------
        ##
	## For a CNSet object (from the crlmm package):
        ##
	##---------------------------------------------------------------------------
	library(oligoClasses)
	library2(crlmm)
	data(cnSetExample, package="crlmm")
	## coerce to an object with log R ratios and B allele frequencies
	oligosetlist <- constructOligoSetListFrom(cnSetExample)
	oligoset <- oligosetlist[[1]]
	res <- hmm(oligoset, p.hom=0, prOutlierBAF=list(initial=1e-4, max=1e-1, maxROH=1e-3))
	res <- res[["NA19007"]]
	rd <- res[state(res)!=3 & numberProbes(res) >= 5, ]
	elementMetadata(rd)$sampleId <- "NA19007"
	if(FALSE){
		## a lattice display for multiple CNV calls ranges.
		library(Biobase)
		library(IRanges)
		xyplotLrrBaf(rd, oligoset,
			     frame=200e3,
			     panel=xypanelBaf,
			     cex=0.5,
			     scales=list(x=list(relation="free"),
			     y=list(alternating=1,
			     at=c(-1, 0, log2(3/2), log2(4/2)),
			     labels=expression(-1, 0, log[2](3/2), log[2](4/2)))),
			     par.strip.text=list(cex=0.7),
			     ylim=c(-3,1),
			     col.hom="grey50",
			     col.het="grey50",
			     col.np="grey20",
			     key=list(text=list(c(expression(log[2]("R ratios")), expression("B allele freqencies")),
				      col=c("grey", "blue")), columns=2))
		frange <- makeFeatureGRanges(oligoset)
		i <- subjectHits(findOverlaps(rd[1,], frange))
		b <- baf(oligoset)[i, 1]
		b <- b/1000
		hist(b, breaks=100)
	}
}
\keyword{methods}

