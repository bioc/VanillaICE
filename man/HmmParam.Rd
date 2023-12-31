% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-HmmParam.R
\name{HmmParam}
\alias{HmmParam}
\alias{HmmParam,missing-method}
\alias{HmmParam,matrix-method}
\alias{show,HmmParam-method}
\alias{nrow,HmmParam-method}
\alias{ncol,HmmParam-method}
\title{Constructor for HmmParam class}
\usage{
HmmParam(
  emission = matrix(0, 0, 0),
  emission_param = EmissionParam(),
  transition = rep(0.99, nrow(emission)),
  chromosome = character(nrow(emission)),
  loglik = LogLik(),
  viterbi = Viterbi(),
  compute_posteriors = TRUE,
  verbose = FALSE
)

\S4method{show}{HmmParam}(object)

\S4method{nrow}{HmmParam}(x)

\S4method{ncol}{HmmParam}(x)
}
\arguments{
\item{emission}{A matrix of emission probabilities}

\item{emission_param}{an object of class \code{EmissionParam}}

\item{transition}{vector of transition probabilities whose length
is N-1, where N is the number of markers.  User should provide the
probability that the state at marker j is the same as the state at
marker j-1.  It is assumed that the probability of transitioning to
state_j from state_j-1 is the same for all states != state_j-1.}

\item{chromosome}{character vector}

\item{loglik}{an object of class \code{LogLik}}

\item{viterbi}{an object of class \code{Viterbi}}

\item{compute_posteriors}{logical}

\item{verbose}{logical}

\item{object}{a \code{HmmParam} object}

\item{x}{a \code{HmmParam} object}
}
\description{
Contains emission probabilities, parameters for emission
probabilities, and transition probabilities required for computing
the most likely state path via the Viterbi algorithm
}
\examples{
HmmParam()
}
