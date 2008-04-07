\name{calls.emission}
\alias{calls.emission}
\title{Calculates emission probabilities for the genotype calls}
\description{
  Emission probabilites for the genotype calls 
}
\usage{
calls.emission(object)
}
\arguments{
  \item{object}{Object inheriting from class \code{HmmOptions}}
}
\details{
  
   If calls.ICE = TRUE in the \code{HmmOptions} object, then the
   emission probabilities for the genotype calls take into account the
   confidence scores for the called genotypes.
   
 }
 
 \value{
   
   Array of emission probabilities for the genotype calls (on the log
   scale)
   
 }
 
\author{R. Scharpf}
\seealso{\code{\link{calls.ICE}}}
\keyword{methods}

