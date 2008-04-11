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


\details{

   If TRUE, the confidence scores (assumed to be inverse standard
      errors) are used by the hidden Markov model.  This has the
      desirable consequence of providing more local smmoothing for
      precise estimates of copy number and more global smoothing for
      copy number estimates with higher standard errors. When FALSE, a
      robust estimate of the chip variability is calculated across all
      of the autosomes.
      
}

\value{
  Logical
}
\author{R. Scharpf}
\seealso{See the VanillaICE vignette}
\keyword{methods}

