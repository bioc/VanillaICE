\name{calculateDistance}
\alias{calculateDistance}
\alias{calculateDistance,SnpLevelSet-method}
\title{Calculates the }
\description{
  Calculates the distance between adjacent SNPs.
}
\usage{
calculateDistance(object)
}
\arguments{
  \item{object}{Object inheriting from \code{SnpLevelSet}}
}
\details{
    Calculates the physical distance between adjacent SNPs along a chromosome.
}
\value{
  A vector of integers
}
\author{R. Scharpf}
\note{ The transition probabilities in the Viterbi algorithm are
  calculated as a function of the physical distance (d) between adjacent
  SNPs.

  probability of staying in state S = exp(-2 * d/100*1e6)
  
}
\seealso{\code{\link{viterbi}}}
\keyword{methods}
\keyword{manip}
