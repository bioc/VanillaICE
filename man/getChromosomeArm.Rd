\name{getChromosomeArm}
\alias{getChromosomeArm}
\title{Get the chromosome arm}
\description{
  Retrieves the chromosome arm (p or q) for each SNP.
}
\usage{
getChromosomeArm(snpset)
}
\arguments{
  \item{snpset}{An object inheriting from class \code{SnpLevelSet}}
}
\value{
  A vector of class 'character'.  
}
\note{

  This function uses data(chromosomeAnnotation) in the SNPchip package
  to retrieve this information.

}
\seealso{\code{\link[SNPchip]{chromosomeAnnotation}}}
\author{RS}
\keyword{manip}

