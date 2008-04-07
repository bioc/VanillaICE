\name{transitionScale}
\alias{transitionScale}
\alias{transitionScale<-}
\title{Scales the transition probabilities of the hidden Markov model}

\description{
  
  Accessor and replacement methods for scaling the transition probabilities
  
}

\usage{
transitionScale(object)
transitionScale(object) <- value
}
\arguments{
  \item{object}{ An object of class \code{HmmParameter}}
  
  \item{value} {A matrix of dimension STATES x STATES, where STATES is
    the number of hidden states.}

}

\details{

  The probability of remaining in the same state, $P(S_t = S_{t+1})$
(the diagonal of the transition probability matrix) is a function of the
distance between SNPs.  The probability of transitioning to some other
state is epsilon, where epsilon = 1 - $P(S_t = S_{t+1})$.  The epsilon
is split among STATES-1 states.  By default, the probability of
transitioning from an altered state back to the normal state is twice as
likely as the probability of transitioning between two altered states.
The weights for epsilon are provided in the \code{transitionScale}
matrix in objects of class \code{HmmParameter}.
  
}

\value{
  A matrix
}

\author{R. Scharpf}
\seealso{ \code{\link{scaleTransitionProbability}}}
\keyword{methods}

