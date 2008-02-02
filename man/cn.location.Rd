\name{cn.location}
\alias{cn.location}
\title{ Accessor for the location of the hidden states for copy number}
\description{
The HMM assumes that conditional on the hidden state, the observed copy
number is approximately Gaussian.  The \code{cn.location} is an accessor
for the mean of the Gaussian distribution for each of the hidden
states.  
}
\usage{
cn.location(object)
}
\arguments{
  \item{object}{ An object of class \code{HmmParameter}}
}
\value{
  A numerical vector.  Typically these values are on the log-scale.
}
\author{R. Scharpf}
\seealso{\code{\link{HmmParameter-class}}}
\examples{
}
\keyword{methods}

