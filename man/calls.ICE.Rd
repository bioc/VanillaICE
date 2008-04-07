\name{calls.ICE}
\alias{calls.ICE}
\alias{calls.ICE<-}
\title{Indicates whether to use confidence scores of the observed copy
  number estimates in the HMM}

\description{ In the current framework, the confidence scores are
  assumed to be the inverse of the standard error }

\usage{
calls.ICE(object)
calls.ICE(object) <- value
}
\arguments{
  \item{object}{An object of class \code{HmmOptions}}
  \item{value}{Logical.  If TRUE, confidence scores of the genotype
    calls are used to calculate emission probabilities}
}
\value{
  Logical
}
\author{R. Scharpf}
\seealso{See the VanillaICE vignette}
\keyword{methods}