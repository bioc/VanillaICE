\name{hmm}
\alias{hmm}
\title{Fits the hidden Markov model}
\description{
Fits the hidden Markov model
}
\usage{
hmm(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HmmParameter}}
  \item{\dots}{Not implemented}
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
params <- new("HmmParameter", snpset=chromosome1,
	      states=c("D", "N", "L", "A"), ##Deletion, Normal, Loss of hetozygosity, Amplification
	      cn.location=c(1, 2, 2, 3),
	      gte.state=c(0.99, 0.8, 0.99, 0.8))
fit <- hmm(params, chromosome1)
}
\keyword{models}
