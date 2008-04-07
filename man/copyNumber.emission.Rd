\name{copyNumber.emission}
\alias{copyNumber.emission}
\title{Calculate emission probabilities for copy number}
\description{
  Calculate emission probabilities for the copy number
}
\usage{
copyNumber.emission(object)
}

\arguments{
  \item{object}{Instance of class \code{HmmOptions}}
}

\details{

     If \code{copyNumber.ICE} = TRUE in the \code{HmmOptions} object,
   then the emission probabilities for the copy number estimates take
   into account the confidence scores stored in the cnConfidence slot of
   the object inherting from class \code{SnpLevelSet}.  The confidence
   scores are assumed to be inverse standard errors.
   
}

\value{

     Array of emission probabilities for the copy number estimates (on
     the log scale)

}
\author{R. Scharpf}
\seealso{\code{\link{copyNumber.ICE}}}
\keyword{methods}

