\name{copyNumberLimits}
\alias{copyNumberLimits}

\title{
  Constraints for updating the means for the copy number states in the
  hidden markov model.
}
\description{
  Constraints for updating the means for the copy number states in the
  hidden markov model.
}
\usage{
copyNumberLimits(is.log)
}
\arguments{
  \item{is.log}{
    \code{logical}: whether the copy number estimates are on the log scale
  }
}
\details{

  Not indented to be called directly -- used by packages that depend on
  VanillaICE.

}

\value{
  A \code{numeric} vector of length 2 giving the lower and upper bounds
  for the copy number estimates.
}

\author{
R. Scharpf
}

\keyword{manip}

