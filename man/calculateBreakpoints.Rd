\name{calculateBreakpoints}
\alias{calculateBreakpoints}
\title{Calculating breakpoints in the hidden states}
\description{
  Computes indices for the breakpoints from the matrix of predicted
  states.
}
\usage{
calculateBreakpoints(x, position, chromosome, states, digits = 2, sampleNames)
}
\arguments{
  \item{x}{vector of predictions}
  \item{position}{Vector of physical positions}
  \item{chromosome}{ Vector of chromosomes (character string)}
  \item{states}{Character string:  the hidden states}
  \item{digits}{ Whether to round the results (e.g., size of altered region)}
  \item{sampleNames}{ Names of the samples}
}
\details{
  None yet
}
\value{
  List of breakpoint matrices.  Each matrix in the list corresponds to
  the matrix of breakpoints identified in a single sample.
}
\author{R. Scharpf}
\seealso{\code{\link{HmmPredict-class}}, \code{\link{hmm}}}
\keyword{manip}

