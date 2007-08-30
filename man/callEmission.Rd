\name{callEmission}
\alias{callEmission}
\title{Calculate emission probabilities for the ICE HMM}

\description{

  When confidence estimates are incorporated in the ICE HMM for
  genotype calls, the observed data at a SNP is the genotype call and
  the uncertainty measure. The joint distribution of the call and the
  uncertainty measure depends on the underlying hidden state.  This
  function calculates the emission probability for the joint sequence
  of calls and confidence scores. See reference for additioal details.

}

\usage{
callEmission(object, P.CHOM.Normal, P.CHOM.LOH, SAMPLE=1)
}

\arguments{
  \item{object}{An object inheriting from \code{\link[SNPchip]{AnnotatedSnpSet}}}
  \item{P.CHOM.Normal}{ probability of called homozygous given the true
  underlying state is normal (retention of heterozygosity)}
  \item{P.CHOM.LOH}{ probability of called homozygous given the true
  underlying state is loss of heterozygosity}
  \item{SAMPLE}{Integer: which sample}
}
\details{

  This function is essentially a wrapper for
  \code{getHapmapProbabilities} that allows the object to contain the
  merged 50k Xba and Hind chips, or the merged 250k Nsp and Sty chips.
  Because the distribution of the uncertainty estimates are
  chip-specific (enzyme-specific), this function computes
  enzyme-specific emission probabilities (note: one of the columns in
  \code{featureData} must have the name 'enzyme' and the elements in
  the vector must be 'Xba' and 'Hind' or 'Nsp' and 'Sty'.).

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
}

\references{ }
\author{}


\note{
}

\seealso{

  \code{\link{getHapmapProbabilities}},   \code{\link[oligo]{crlmm}},
    \code{\link[Biobase]{featureData}}
  For the callsConfidence package, see
  \url{http://www.biostat.jhsph.edu/~rscharpf/software/index.html}

}

\examples{
  data(chromosome1)
  ##Create an object of class AnnotatedSnpCallSet
  callset <- new("AnnotatedSnpSet",
                 calls=calls(chromosome1),
                 callsConfidence=callsConfidence(chromosome1),
                 featureData=featureData(chromosome1),
                 phenoData=phenoData(chromosome1),
                 annotation=annotation(chromosome1))
  ##1 = homozygous (AA or BB), 2 = heterozygous
  calls(callset)[calls(callset) == 3] <-  1
  ##make sure that enzyme information is stored in featureData

\dontrun{
  require("callsConfidence")
    emissionProb <- callEmission(object=callset,
                                 P.CHOM.Normal=0.70,
                                 P.CHOM.LOH=0.999)
    range(emissionProb)
}
}



\keyword{manip}

