\name{genotypeEmissionCrlmm}
\alias{genotypeEmissionCrlmm}
\title{Estimate the emission probabilities using confidence scores from CRLMM}
\description{
  Estimate the emission probabilities that incorporate information on
  the confidence scores for the genotype calls. 
}
\usage{
genotypeEmissionCrlmm(genotypes, conf, pHetCalledHom = 0.001, pHetCalledHet = 0.995, pHomInNormal = 0.99, pHomInRoh = 0.999, annotation)
}
\arguments{
  \item{genotypes}{Matrix of genotypes}
  \item{conf}{Matirx of confidence scores (see details).}
  \item{pHetCalledHom}{The probability that a truly heterozygous SNP is
    incorrectly called homozygous (incorrect call).}
  \item{pHetCalledHet}{The probability that a truly heterozygous SNP is
    called heterozygous (correct call).}
  \item{pHomInNormal}{ The probability of a homozygous genotype call in
    the 'normal' state.}
  \item{pHomInRoh}{The probability of a homozygous genotype call in a
    region of homozygosity.}
  \item{annotation}{The cdf name (e.g., "genomewidesnp6")}
}
\details{
  The confidence scores by crlmm are saved as an integer: 1000*log(1-p),
  where p is the probability that the genotype call is correct.

  The reference distribution of confidence scores are available for the following
  Affymetrix platforms: affy6, nsp250, and sty250k.
  
}

\value{
  An R x C x X array of emission probabilities, where

  R = number of loci (SNPs)
  C = number of samples
  S = number of states
  
}
\references{RB Scharpf et al. (2008), Annals of Applied Statistics}
\author{R Scharpf}
\keyword{htest}

