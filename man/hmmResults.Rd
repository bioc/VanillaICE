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

\seealso{\code{\link{xyplot}}, \code{\link{hmm}}}

\examples{
	data(hmmResults)
	class(hmmResults)
}
\keyword{datasets}
