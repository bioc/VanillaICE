\name{oligoSetList-methods}
\docType{methods}
%\alias{hmm}
%\alias{hmm-methods}
\alias{[[,oligoSetList-method}
\alias{dims,oligoSetList-method}

\title{Methods for oligoSetList class}

\description{

  The oligoSetList class is a container for genotypes, B allele
  frequencies, and copy number organized by chromosome.  Each element in
  the list class contains low-level summaries and phenotypic information
  for a single chromosome. The organization by chromosome facilitates
  parallelization of methods to identify copy number alterations.

}

\section{Methods}{

  For each of the following methods, \code{object} is an instance of
  class \code{oligoSetList}.


  \describe{

    \item{}{\code{object[[i]]}:

      \code{i} must be an integer. Return a \code{oligoSnpSet} object
      for the ith element in the \code{oligoSetList} object.

    }

    \item{}{\code{object[i]}:

      \code{i} can be a vector of integers. Returns an object of the same
      class with length equal to the length of the \code{i} vector.
    }

    \item{}{dims(object):

      Return object dimensions

    }
  }
}

\seealso{\code{\link{hmmBeadStudioSetList}}, \code{\link{hmmOligoSetList}}}

\examples{
	library(oligoClasses)
	library2(crlmm)
	data(cnSetExample, package="crlmm")
	## coerce to an object with log R ratios and B allele frequencies
	oligosetlist <- constructOligoSetListFrom(cnSetExample)
	oligoset <- oligosetlist[[1]]
}
\keyword{methods}

