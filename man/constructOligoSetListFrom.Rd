\name{constructOligoSetListFrom}
\alias{constructOligoSetListFrom}
\alias{constructBafLrrSetListFrom}

\title{Construct a chromosome-list container from a CNSet object}
\description{
  The oligoSetList contains genotype calls, genotype call probabilities,
  BAFs, and log R ratios.
}
\usage{
constructOligoSetListFrom(object, ...)
constructBafLrrSetListFrom(object, ...)
}
\arguments{
  \item{object}{
    A \code{CNSet} object.
}
  \item{\dots}{
    Additional arguments to \code{calculateRBaf}.  In particular,
    \code{batch.name} (name of batch) and \code{chrom} (which
    chromosomes). If \code{batch.name} and \code{chrom} are missing, a \code{oligoSetList}
    containing all batches and chromosomes 1-22 and X are processed.
}
}

\value{
An object of class \code{oligoSetList}
}
\author{
R. Scharpf
}

\seealso{
\code{\link{calculateRBaf}}
}
\examples{
library(oligoClasses)
library2(Biobase)
data(cnSetExample, package="crlmm")
constructOligoSetListFrom(cnSetExample)
constructBafLrrSetListFrom(cnSetExample)
}
\keyword{manip}

