\name{hmm}
\alias{hmm}
\title{Fits the hidden Markov model}
\description{
Fits the hidden Markov model
}
\usage{
hmm(object, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HmmOptions}}
  \item{params}{An object of class \code{HmmParameter}}
}
\details{
  None yet.
}
\value{
  An object of class \code{HmmPredict}
}

\references{

    RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

  }
\author{R. Scharpf}
\seealso{

  \code{\link{HmmParameter-class}}
  \code{\link{HmmPredict-class}}
  \code{\link[oligoClasses]{SnpLevelSet-class}}

}
\examples{
data(chromosome1)
chromosome1 <- chromosome1[1:500, ]
options <- new("HmmOptions",
	       states=c("D", "N", "L", "A"),
	       snpset=chromosome1,
	       copyNumber.location=c(1, 2, 2, 3),
	       probHomCall=c(0.99, 0.7, 0.99, 0.7))
validObject(options)
params <- new("HmmParameter", 
	      states=states(options),
	      initialStateProbability=0.99)
cn.emission <- copyNumber.emission(options)
gt.emission <- calls.emission(options)
emission(params) <- cn.emission + gt.emission ##log scale
genomicDistance(params) <- exp(-2 * calculateDistance(options)/(100*1e6))
transitionScale(params) <- scaleTransitionProbability(options)
validObject(params)
fit <- hmm(options, params)
}
\keyword{models}
