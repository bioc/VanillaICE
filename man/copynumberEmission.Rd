\name{copynumberEmission}
\alias{copynumberEmission}
\title{Emission probabilities for copy number}
\description{
Emission probabilities for copy number
}
\usage{
copynumberEmission(copynumber, states, mu, uncertainty, takeLog, verbose = TRUE)
}
\arguments{
  \item{copynumber}{class matrix}
  \item{states}{class character}
  \item{mu}{class numeric: mean of hidden states for Gaussian}
  \item{uncertainty}{uncertainty of copy number estimates}
  \item{takeLog}{whether to take the log of the copy number AND the mu argument}
  \item{verbose}{logical}
}
\details{
}
\value{
  \item{array}{Array of emission probabilities.  Dimension 1: SNPs, Dimension 2: samples, Dimension3: states}
}
\seealso{}
\examples{
}
\keyword{methods}

