\name{scaleTransitionProbability}
\alias{scaleTransitionProbability}
\title{Scales the probability of transitioning between two states in the
  Viterbi algorithm}

\description{
Function only returns a matrix of 1's (by default, no scaling of the
transition probabilities).  This may be changed in future releases.
}
\usage{
scaleTransitionProbability(states, SCALE = 2, normalLabel = "N")
}
\arguments{
  \item{states}{character string of hidden states}
  \item{SCALE}{Numeric. See details.}
  \item{normalLabel}{Character: label used for the baseline state.}
}
\details{
}
\value{
  A matrix of dimension S x S, where S is the number of states.  The
  diagonal should always be 1 (we do not scale the probability of
  remaining in the same state).  A slot for storing this matrix is
  provided in the \code{HmmParameter} class.  
}
\author{RS}
\seealso{\code{\link{genomicDistance}}, \code{\link{transitionScale}}, \code{\link{HmmParameter-class}}}
\examples{
}
\keyword{manip}

