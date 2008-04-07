\name{genomicDistance}
\alias{genomicDistance}
\alias{tau}
\alias{genomicDistance<-}
\title{A rough estimate of the genomic distance between SNPs}

\description{
  Calculates the genomic distance between SNPs as a function of the
  physical distance.  This distance is used for calculating transition
  probabilities in the viterbi algorithm.
}

\usage{
genomicDistance(object)
tau(object)
genomicDistance(object) <- value
}

\arguments{
  \item{object}{ An object inheriting from \code{HmmOptions}}
  \item{value}{A vector of length T - 1, where T is the number of SNPs}
}
\details{

  The genomic distance is calculated as exp(-2*d), where d is the
  physical distance in 100MB units.

  In the Viterbi algorithm, the probability of remaining in state S at
  SNP t+1 is the genomic distance.  The probability of transitioning to
  one of the alternative states is 1 - the genomic distance.  The
  \code{transitionScale} function can be used to give more weight to
  transitioning back to the baseline ('normal') state.
  
}

\value{
    
    A numeric vector with length T-1.

}

\author{R. Scharpf}
\seealso{\code{\link{transitionScale}}}
\keyword{methods}