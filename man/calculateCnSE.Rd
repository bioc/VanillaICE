\name{calculateCnSE}
\alias{calculateCnSE}
\title{Calculates copy number standard errors}
\description{
Calculates SNP-specific standard errors for copy number using a reference set.
}
\usage{
calculateCnSE(object, referenceSet, epsilon = 0.1)
}
\arguments{
  \item{object}{ An object of class \code{SnpCopyNumberSet} or
    \code{oligoSnpSet}.  }
  \item{referenceSet}{ An object of class \code{SnpCopyNumberSet} or
    \code{oligoSnpSet}. If missing, this function uses the samples in
    the \code{object} to estimate across sample standard deviations }
  \item{epsilon}{minimum standard error}
}
\details{

  Calculates SNP-specific standard errors from a reference distribution
  as the product of the within copy number standard deviation of the
  sample (test set) and the across samples standard deviation from a
  reference set.  The across samples standard deviation is scaled by the
  median (e.g., across sample sd = across sample sd/ median(across
  sample sd)).
  
}
\value{
  A matrix of standard errors
}
\author{R. Scharpf}
\keyword{manip}

