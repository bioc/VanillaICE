\name{hapmapProbabilities}
\alias{hapmapProbabilities}
\title{CRLMM confidence scores on HapMap data}
\description{
  
  Wrapper for extracting the chip-specific distribution of CRLMM
  confidence scores for HapMap data where we assume we know the true
  genotype. Presently only supported for Affymetrix 50k and 250k SNP
  chips and onnly relevant for fitting the ICE HMM (see reference
  below).
  
}

\usage{
hapmapProbabilities(x)
}

\arguments{
  \item{x}{Name of chip. Must be 'pd.mapping...'}
}

\details{

  Confidence scores in a \code{SnpCallSet} created by \code{crlmm} are
  stored in the element "pAcc" in the assay data.

}

\value{
  Matrix of confidence scores
}

\references{

  RB Scharpf et al. (2007), Department of Biostatistics
  Working Papers, Working Paper 136

}

\author{ }
\note{ 
}

\seealso{

  \code{\link[oligo]{crlmm}},
  \code{\link[oligo]{SnpCallSet-class}},
  \code{\link[Biobase]{AssayDataElement}}

}

\examples{
}

\keyword{data}

