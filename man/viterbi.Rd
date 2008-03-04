\name{viterbi}
\alias{viterbi}
\title{viterbi algorithm}

\description{

  The Viterbi algorithm for computing the most likely state sequence
  given a model

}
\usage{
viterbi(beta, pi, tau, states, arm, tau.scale)
}

\arguments{
  \item{beta}{emission probabilities (log scale)}
  \item{pi}{initial state probabilities (log scale)}
  \item{tau}{transition probabilities (log scale)}
  \item{states}{character string of hidden states}
  \item{arm}{character string of chromosomal arm}
  \item{tau.scale}{scales the probability of transitioning between an
  altered state to a normal state relative to the probability of
  transitioning between two altered states.}
}
\details{

  The Viterbi algorithm is fit independently to each chromosomal arm.

  Argument \code{tau.scale} is a matrix that scales the probability of
  transitioning from an altered state to a normal state to the
  probability of transitioning between two altered states.

  When initializing an object of class \code{HmmParameter}, the default
  is to make the probability of transitioning from an altered state to a
  normal state twice the probability of transitioning between two
  altered states. Further flexibility could be obtained by replacing the
  tau.scale matrix in an instance of the \code{HmmParameter} class with
  your own matrix.  The dimension of this matrix must be S x S, where S
  is the number of hidden states.

}
\value{
  vector of integers
}
\author{R. Scharpf}
\keyword{models}

