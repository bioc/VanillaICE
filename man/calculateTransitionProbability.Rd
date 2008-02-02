\name{calculateTransitionProbability}
\alias{calculateTransitionProbability}
\title{Calculates transition probabilities}
\description{
Calculates transition probabilities as a function of the distance
between SNPs
}
\usage{
calculateTransitionProbability(object, options, ...)
}

\arguments{
  \item{object}{ An object inheriting from \code{SnpLevelSet} (e.g., an \code{oligoSnpSet})}
  \item{options}{An object of class \code{HmmOptions}}
  \item{\dots}{Not implemented }
}
\details{

  The transition probabilities are SNP-specific.  We assume that the
  probability of remaining in the same state at SNP t+1 is 

  exp(-2*d)

  d = the distance in 100 MB.

  The probability of transitioning to a different state at SNP t+1 is
  divided equally among the other states.  To give more weight to
  transitioning between certain states (e.g., for going back to the
  normal state or prohibiting transitions between two altered states),
  see \code{tau.scale}.
  }

  \value{
    
A matrix with S*S columns (S = number of hidden states) and R rows (R =
number of SNPs)

}

\author{R. Scharpf}
\seealso{\code{\link{tau.scale}}}
\keyword{methods}
