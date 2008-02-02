\name{arm}
\alias{arm}
\title{Accessor for chromosome arm}
\description{
  Accessor for chromosome arm
}
\usage{
arm(object)
}
\arguments{
  \item{object}{ an object of class \code{HmmParameter}}
}

\value{
  A character string.
}
\author{R. Scharpf}
\seealso{\code{\link{HmmParameter-class}}}
\examples{
data(chromosome1)
chromosome1 <- chromosome1[1:10,]
params <- new("HmmParameter", snpset=chromosome1)
arm(params)[1:5]
}
\keyword{methods}

