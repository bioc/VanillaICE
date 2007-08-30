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
  \item{\dots}{additional arguments to the HMM}
}
\details{

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
  fitVanilla <- hmm(hmmset, flavor="vanilla", arm="both")

  require(callsConfidence))
    fitIce <- hmm(hmmset, flavor="ICE", arm="both")
}
}
\keyword{models}
\keyword{smooth}
