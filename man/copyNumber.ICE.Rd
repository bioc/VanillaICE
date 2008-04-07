\name{copyNumber.ICE}
\alias{copyNumber.ICE}
\alias{copyNumber.ICE<-}
\title{Indicates whether to use confidence scores of the observed copy
  number estimates in the HMM}

\description{ In the current framework, the confidence scores are
  assumed to be the inverse of the standard error }

\usage{
copyNumber.ICE(object)
copyNumber.ICE(object) <- value
}
\arguments{
  \item{object}{An object of class \code{HmmOptions}}
  \item{value}{Logical. If TRUE, confidence scores of the copy number
    estimates are used to calculate emission probabilities.}
}
\value{
  Logical
}
\author{R. Scharpf}
\seealso{See the VanillaICE vignette}
\keyword{methods}

