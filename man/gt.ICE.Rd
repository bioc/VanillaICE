\name{gt.ICE}
\alias{gt.ICE}
\alias{gt.ICE<-}
\title{Indicates whether to use confidence scores of the observed copy
  number estimates in the HMM}

\description{ In the current framework, the confidence scores are
  assumed to be the inverse of the standard error }

\usage{
gt.ICE(object)
}
\arguments{
  \item{object}{An object of class \code{HmmParameter}}
}
\value{
  Logical
}
\author{R. Scharpf}
\seealso{See the VanillaICE vignette}
\keyword{methods}