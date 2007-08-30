\name{pCHOM}
\alias{pCHOM}
\alias{pCHOM<-}
\title{Probability of a homozygous genotype calls}
\description{
  
  Accessor for the probability of a homozygous genotype calls
  conditional on the hidden state
  
}
\usage{
pCHOM(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{HmmSnpCallSet} or \code{HmmSnpSet}}
}
\details{
}
\value{
  
  A numerical vector of probabilities with length equal to the number of
  hidden states
  
}
\references{}
\author{}
\note{ 
}
\seealso{}
\examples{
data(chromosome1)
pCHOM(chromosome1)
}
\keyword{methods}

