\name{getTransitionProbability}
\alias{getTransitionProbability}
\title{Transition probabilities for the HMM}

\description{
  
  The default transition probabilities for the HMM can be a matrix with
  dimension S x S, where S is the number of hidden states.
  Alternatively, one may specify transition probabilities that are a
  function of the distance (in base pairs) between adjacent SNPs.
  
}

\usage{
getTransitionProbability(object)
}
\arguments{
  \item{object}{ An object of one of the HMM classes defined in \code{ICE}}
}
\details{
}
\value{

  A transition probability matrix.
  
}

\references{RB Scharpf et al. (2007), Department of Biostatistics
  Working Papers, Working Paper 136 and CITATION FOR MA NEWTON PAPER}

\author{}
\note{ 
}
\seealso{\code{\link{distance}}}
\examples{
}
\keyword{attribute}

