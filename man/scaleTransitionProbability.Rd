\name{scaleTransitionProbability}
\alias{scaleTransitionProbability}
\title{Scales the probability of transitioning between two states in the
  Viterbi algorithm}

\description{
Scales the probability of transitioning between two states in the
  Viterbi algorith
}
\usage{
scaleTransitionProbability(object, SCALE = 2, normalLabel = "N")
}
\arguments{
  \item{object}{Object of class \code{HmmOptions}}
  \item{SCALE}{Numeric. See details.}
  \item{normalLabel}{Character: label used for the baseline state.}
}
\details{
  If SCALE is 2, the probability of transitioning from an altered state
  back to the normal state is twice the probability of transition from
  an altered state to another altered state.
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
data(chromosome1)
options <- new("HmmOptions",
	       states=c("D", "N", "L", "A"),
	       snpset=chromosome1,
	       copyNumber.location=c(1, 2, 2, 3),
	       probHomCall=c(0.99, 0.7, 0.99, 0.7))
scaleTransitionProbability(options)
}
\keyword{manip}

