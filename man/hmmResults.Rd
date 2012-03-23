\name{hmmResults}
\alias{hmmResults}
\docType{data}
\title{
  Example output from hmm
}
\description{
  Example output from \code{hmm} method applied to simulated data.
}
\usage{data(hmmResults)}
\format{
  A \code{RangedDataHMM} object.
}
\details{
  The results of a 6-state HMM fit to simulated copy number and genotype
  data.
}

\seealso{\code{\link{xyplot}}}

\examples{
	library(oligoClasses)
	data(oligoSetExample, package="oligoClasses")
	oligoSet <- oligoSet[chromosome(oligoSet) == 1, ]
	hmmResults <- hmm(oligoSet)
	state(hmmResults)
	##
	## Plotting ranges:
	##
	if(require(SNPchip) && require(IRanges)){
		## Plot the data for the second range with a blue
		## border, and frame the region by 10 Mb on each side
		## of the state boundary.
		##
		xyplot(cn~x, oligoSet, range=hmmResults[2, ], frame=10e6,
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
		xyplot(cn~x | range, oligoSet, range=hmmResults, frame=2e6,
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
}
\keyword{datasets}
