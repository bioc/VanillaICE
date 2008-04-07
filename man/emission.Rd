\name{emission}
\alias{emission}
\alias{emission<-}
\title{Accessor for emission probabilities}
\description{
  Accessor for emission probabilities for the hidden Markov model
}
\usage{
emission(object)
emission(object) <- value
}
\arguments{
  \item{object}{ An object of class \code{HmmParameter}}
  \item{value}{An array of emission probabilities (log scale)}
}

\value{
  An array.
}

\author{R. Scharpf}
\seealso{\code{\link{HmmParameter-class}}, \code{\link{HmmOptions-class}}}
\keyword{methods}

