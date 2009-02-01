\name{copynumberEmission}
\alias{copynumberEmission}
\title{Emission probabilities for copy number}
\description{
Emission probabilities for copy number
}
\usage{
copynumberEmission(copynumber, states, mu, sds, takeLog, verbose = TRUE)
}
\arguments{
  \item{copynumber}{matrix}
  \item{states}{character string}
  \item{mu}{numeric: mean of hidden states for Gaussian}
  \item{sds}{standard deviations of copy number estimates}
  \item{takeLog}{logical: if TRUE, this function takes the log of the
    copy number AND mu arguments to this function}
  \item{verbose}{logical}
}
\details{
}
\value{
  
  \item{array}{Array of emission probabilities on the log scale.
    Dimension 1: SNPs, Dimension 2: samples, Dimension3: states}
  
}
\seealso{}
\examples{
}
\keyword{methods}

