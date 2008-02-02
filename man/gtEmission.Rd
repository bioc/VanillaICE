\name{gtEmission}
\alias{gtEmission}
\title{Calculates emission probabilities for the genotype calls}
\description{
  Emission probabilites for the genotype calls 
}
\usage{
gtEmission(x, options, ...)
}
\arguments{
  \item{x}{Vector of integers for the observed genotypes }
  \item{options}{ An object of class \code{HmmOptions}}
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
   if gt.ICE = TRUE in the \code{HmmOptions} object, then the emission
   probabilities for the genotype calls take into account the confidence
   scores for the observed genotypes
}
\value{
  Vector of the logged emission probabilities
}
\author{R. Scharpf}
\seealso{\code{\link{cn.ICE}}}
\keyword{methods}

