\name{copynumberEmission}
\alias{copynumberEmission}
\title{Emission probabilities for copy number}
\description{
Emission probabilities for copy number
}
\usage{
copynumberEmission(copynumber, states, mu, sds, takeLog, verbose = TRUE,
na.rm=TRUE)
}
\arguments{
  \item{copynumber}{matrix}
  \item{states}{character string}
  \item{mu}{numeric: mean of hidden states for Gaussian}
  \item{sds}{standard deviations of copy number estimates}
  \item{takeLog}{logical: if TRUE, this function takes the log of the
    copy number AND mu arguments to this function}
  \item{verbose}{logical}
  \item{na.rm}{The default is to ignore missing values when calculating robust standard deviations}
}
\value{
  
  \item{array}{Array of emission probabilities on the log scale.
    Dimension 1: SNPs, Dimension 2: samples, Dimension3: states}
  
}
\details{

  By default, this func estimates the scale parameter for the Normal
  distribution from the supplied data using the median absolute
  deviation (MAD). However, different standard deviations can be
  supplied by the user with the argument \code{sds}. The supplied
  standard deviations must be of the same dimension as the copy number
  matrix. 

}
\seealso{\code{\link{genotypeEmission}}, \code{\link{genotypeEmissionCrlmm}}}
\keyword{methods}

