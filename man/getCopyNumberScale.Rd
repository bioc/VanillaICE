\name{getCopyNumberScale}
\alias{getCopyNumberScale}
\title{Scale parameter used for Gaussian distribution of copy number}
\description{
  We typically assume that after appropriate transformation the
  distribution of copy number is approximately Gaussian.  This function
  calculates a robust estimate of the spread based on quantiles of the
  observed copy number distribution across all autosomes.
}
\usage{
getCopyNumberScale(object, ...)
}
\arguments{
  \item{object}{Object of class \code{HmmSnpSet} or
    \code{HmmSnpCopyNumberSet}}
  \item{\dots}{additional arguments for subsetting the samples}
}
\details{
}
\value{
  A vector of scales with length equal to the number of hidden states
  (typically 3 or 4)
}

\references{RB Scharpf et al. (2007), Department of Biostatistics
  Working Papers, Working Paper 136}

\author{}
\note{ 
}
\seealso{}
\examples{
}
\keyword{attribute}

