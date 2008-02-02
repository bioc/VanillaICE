\name{tau.scale}
\alias{tau.scale}
\title{Scales the transition probabilities of the hidden Markov model}
\description{
  Scales the transition probabilities of the hidden Markov model
}
\usage{
tau.scale(object)
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
\seealso{ \code{\link{calculateTransitionProbability}}}
\examples{

data(chromosome1)
chromosome1 <- chromosome1[1:50, ]
params <- new("HmmParameter", snpset=chromosome1)

##By default we are twice as likely to transition to a normal state
##than between two altered states
tau.scale(params)

##To make it 10 times as likely to transition from an altered state to
##a normal state
params <- new("HmmParameter", snpset=chromosome1, SCALE=10)

}
\keyword{methods}

