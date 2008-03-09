\name{calculateBreakpoints}
\alias{calculateBreakpoints}
\title{ Calculate breakpoints for altered states from the prediction
  matrix}
\description{
 Calculate breakpoints for altered states from the prediction matrix in
 an object of class \code{HmmPredict}
}
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
\seealso{\code{\link{summary}}}
\keyword{manip}
\keyword{methods}
