\name{getCopyNumberLocation}
\alias{getCopyNumberLocation}
\title{Typically the mean copy number for the hidden states}

\description{

  Use to provide default values of the copy number for the hidden
  states.  We typically assume that the distribution of the observed
  copy number (after suitable transformation) is approximately Gaussian
  and centered at the true copy number.
  
}

\usage{
getCopyNumberLocation(object)
}
\arguments{
  \item{object}{An object of class \code{HmmSnpSet} or \code{HmmSnpCopyNumberSet}}
}
\details{

    By default, we assume that the true underlying copy number for a
  deletion is 1, 2 for a normal region (diploid), and 3 for an
  amplification.

}

\value{
  numerical vector of locations for the Gaussian distribution for copy number
}

\references{}
\author{ who you are }
\note{ Not typically called directly
}
\seealso{ }
\examples{
}
\keyword{attribute}

