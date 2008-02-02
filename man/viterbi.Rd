\name{viterbi}
\alias{viterbi}
\title{viterbi algorithm}

\description{

  The Viterbi algorithm for computing the most likely state sequence
  given a model

}
\usage{
viterbi(beta, pi, tau, states, arm)
}

\arguments{
  \item{beta}{emission probabilities (log scale)}
  \item{pi}{initial state probabilities (log scale)}
  \item{tau}{transition probabilities (log scale)}
  \item{states}{character string of hidden states}
  \item{arm}{character string of chromosomal arm}
}
\details{

  The viterbi algorithm is fit independently to each chromosomal arm.

}
\value{
  vector of integers
}
\author{R. Scharpf}
\keyword{models}

