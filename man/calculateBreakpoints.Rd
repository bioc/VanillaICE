\name{calculateBreakpoints}
\alias{calculateBreakpoints}
\title{ Calculate breakpoints for altered states from the prediction
  matrix}
\description{Wrapper for findBreaks}
\usage{
calculateBreakpoints(object, ...)
}
\arguments{
  \item{object}{Object of class \code{HmmPredict}}
  \item{\dots}{ Not implemented}
}
\details{
  Returns a matrix with information on size of region (Mb), predicted state,
  number of SNPs in region, chromosome number, and position.
}
\value{
  A matrix
}
\author{R. Scharpf}
\seealso{\code{\link{findBreaks}}}
\keyword{manip}
\keyword{methods}
