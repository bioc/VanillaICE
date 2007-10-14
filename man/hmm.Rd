\name{hmm}
\alias{hmm}
\title{Fits hidden Markov models to high throughput SNP chip data.}

\description{

Fits a hidden Markov model to high throughput SNP chip data

}

\usage{
hmm(object, ...)
}
\arguments{
  \item{object}{ an object extending \code{hSet}}
  \item{\dots}{See details}
}
\details{

  \item verbose:  if TRUE, print information about the processing
  
  \item log2: if TRUE, a log2 transformation of copy number estimates
  is performed prior to fitting the hmm.

  \item arm: which arm of the chromosome to fit the hidden markov
  model (p, q, or both). Default is to fit to both arms.

  \item MISSING.CODE: missing code for genotype calls (default is NA)

  \item keepEmission: return emission probabilities  (not recommended)

  \item keepLikelihod: return likelihood (not recommended)

  \item keepTransitionProb: return transition probabilities (not recommended)

  \item SCALE=1: whether to scale the transition probabilities


 }
 
 \value{
   
  An object of the same class
  
}

\references{

  \cite{RB Scharpf et al. (2007), Department of Biostatistics
    Working Papers, Working Paper 136}

}

\author{}

\note{ further notes }

\seealso{

  The package \pkg{callsConfidence} is available here:
  \url{http://www.biostat.jhsph.edu/~rscharpf/software/index.html}

}


\examples{

data(chromosome1)

\dontrun{
  callsIce(hmmset) <- FALSE
  copyNumberIce(hmmset) <- FALSE
  fitVanilla <- hmm(hmmset)
  require(callsConfidence)
  callsIce(hmmset) <- TRUE
  copyNumberIce(hmmset) <- TRUE
  fitIce <- hmm(hmmset, flavor="ICE", arm="both")
}
}
\keyword{models}
\keyword{smooth}
