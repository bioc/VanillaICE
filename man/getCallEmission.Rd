\name{getCallEmission}
\alias{getCallEmission}

\title{Calculates emission probabilities for genotype calls}

\description{
  
Calculates emission probabilities for genotype calls. Only used by the
ICE HMM.

}

\usage{

getCallEmission(object, CALL.STATES = 2, SAMPLE = 1, P.HOM.CHOM = 0.999,
P.HET.CHOM = 1 - 0.999, P.HET.CHET = 0.999, P.HOM.CHET = 1 - 0.999,
P.CHOM.Normal, P.CHOM.LOH, hapmapP)

}

\arguments{
  \item{object}{ Typically an object of class \code{HmmSnpCallSet} }
  \item{CALL.STATES}{ Integer: number of hidden states (default is 2
  -- loss of heterozygosity and retention of heterozygosity)}
  \item{SAMPLE}{Integer: which sample to extract genotype calls from the
  \code{HmmSnpCallSet} object}
  \item{P.HOM.CHOM}{ Numeric. Probability that the true genotype is homozygous
  (AA or BB)
  given that the observed genotype was homozygous (AA or BB)}
  \item{P.HET.CHOM}{ Numeric. Probability that the true genotype is
  heterozygous (AB) given that the observed genotype was homozygous.}
  \item{P.HET.CHET}{ Numeric. Probability that the true genotype is
  heterozygous (AB) given that the observed genotype was heterozygous. }
  \item{P.HOM.CHET}{ Numeric. Probability that the true genotype is
  homozygous given that the observed genotype was heterozygous. }
  \item{P.CHOM.Normal}{ Numeric. Probability that the observed genotype
  is homozygous if the true hidden state is retention of heterozygosity.}
  \item{P.CHOM.LOH}{Numeric. Probability that the observed genotype is
  homozygous given that the true hidden state is loss of heterozygosity.}
  \item{hapmapP}{ A matrix estimated from 269 HapMap samples.  See details.}
}

\details{

  We group the CRLMM probabilities of a correct genotype call for the
  269 HapMap samples into 1 of 4 possible categories:
  
  1:  genotype call is HET, truth is HOM
  2:  genotype call is HET, truth is HET
  3:  genotype call is HOM, truth is HOM
  4:  genotype call is HOM, truth is HET
  

}

\value{

  Returns matrix of emission probabilities that are a function of the
  confidence scores provided by \code{crlmm} for genotype calls.

}

\references{

  \cite{RB Scharpf et al. (2007), Department of Biostatistics Working Papers,
  Working Paper 136 }

}

\author{}

\note{

  Not typically called directly.

  
}

\seealso{

  The R package \pkg{callsConfidence}:
  \url{http://www.biostat.jhsph.edu/~rscharpf/software/index.html}
  
  \code{\link{callEmission}}, \code{\link{hmm}},
  \code{\link[oligo]{crlmm}}

}

\examples{
}

\keyword{manip}

