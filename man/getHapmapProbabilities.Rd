\name{getHapmapProbabilities}
\alias{getHapmapProbabilities}
\title{Estimate probability of correct genotype call}

\description{

This function is used by the ICE HMM described in the reference below
and currently only supported for 50k and 250k Affymetrix SNP chips.  A
nonparametric estimate of the probability that the genotype call is
correct is based on the distribution of \code{crlmm} confidence scores
in the reference HapMap population, where we assume that we know the
true genotype.  A similar approach could be used to estimate
probabilities with other genotyping algorithms.  See reference below for
more information.
  
 }

\usage{
getHapmapProbabilities(object)
}

\arguments{
  \item{object}{ One of the HMM classes defined in \code{ICE}}
}

\details{

  The annotation slot in the object must contain one of the following
  strings that identifies the Affymetrix chip:

  mapping100k (if the 50k xba and 50k hind results have been merged
  into one object)

  mapping50k.xba240  (50k Xba platform)

  mapping50k.hind240 (50k Hind)

  mapping500k  (if the two 250k chips -- Nsp and Sty -- have been merged)

  mapping250k.nsp

  mapping250k.sty

 }
 
 \value{
   
   A vector of probabilities that the genotype call is correct.
   
}


\references{

\cite{  RB Scharpf et al. (2007), Department of Biostatistics
  Working Papers, Working Paper 136}

}

\author{}

\note{

  Estimated probabilities are chip-specific.  This function is typically
  not called directly. 
  
}

\seealso{

  \code{\link{callEmission}}, \code{\link[oligo]{crlmm}}

}
  
\examples{
}

\keyword{nonparametric}

