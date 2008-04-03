\name{transitionScale}
%\alias{tau.scale}
\alias{transitionScale}
\alias{transitionScale<-}
\title{Scales the transition probabilities of the hidden Markov model}
\description{
  Scales the transition probabilities of the hidden Markov model
}
\usage{
%tau.scale(object)
transitionScale(object)
}
\arguments{
  \item{object}{ An object of class \code{HmmParameter}}
}
\details{

  The probability of remaining in the same state, $P(S_t = S_{t+1})$ (the
diagonal of the transition probability matrix) is a function of the
distance between SNPs.  The probability of transitioning to some other
state is epsilon, where epsilon = 1 - $P(S_t = S_{t+1})$.  The epsilon is
split among STATES-1 states (where STATES is the number of states).  By
default, the probability of transitioning from an altered state back to
the normal state is twice as likely as the probability of transitioning
between two altered states.  The weights for epsilon are provided in the
\code{tau.scale} matrix in objects of class \code{HmmParameter}
  
}
\value{
  A matrix
}
\author{R. Scharpf}
%\seealso{ \code{\link{calculateTransitionProbability}}}
\keyword{methods}

