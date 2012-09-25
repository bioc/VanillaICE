\name{SetList-methods}
\Rdversion{1.1}
\docType{methods}
\alias{[,BeadStudioSetList-method}
\alias{[,gSetList-method}
%\alias{[[,gSetList-method}
\alias{[[,BeadStudioSetList,ANY,ANY-method}
\alias{[[,BafLrrSetList,ANY,ANY-method}
\alias{[[<-,BafLrrSetList,ANY,ANY,BafLrrSet-method}
\alias{[[<-,BafLrrSetList,ANY,ANY,BafLrrSet}
\alias{[[,oligoSetList,ANY,ANY-method}
\alias{$,gSetList-method}
\alias{$<-,gSetList-method}
\alias{length,gSetList-method}
\alias{dims,gSetList-method}

\title{BeadStudioSetList methods}

\description{
  Methods for BeadStudioSetList objects
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("BeadStudioSetList", assayDataList, logRRatio, BAF, featureDataList, chromosome, phenoData, annotation, genomeBuild, ...)}.
}

\section{Methods}{

  For the following methods, \code{object} can be a
  \code{BeadStudioSetList} or \code{oligoSetList} instance.

  \describe{

    \item{}{\code{object[i]}:

      Returns an object of the same class as \code{object} with length
      equal to \code{length(i)}.

    }

    \item{}{\code{object[[i]]}:

      Returns a \code{BeadStudioSet} or a \code{oligoSnpSet} object,
      depending on whether the class of \code{object} is a \code{BeadStudioSetList}
      or an \code{oligoSetList}.

    }

    \item{}{\code{object[[i]] <- value }:

      Replaces the ith element of the \code{BafLrrSetList} object by
      \code{value}. The object \code{value} must be a \code{BafLrrSet} object.

    }

    \item{}{\code{object$NAME}, \code{object$NAME <- value}:

      Get or set values for for column NAME in \code{phenoData}.  For
      the get method, \code{NAME} must be an element of
      \code{varLabels(object)}.  \code{value} must be the same length as
      \code{ncol(object)}.

    }

    \item{}{\code{hmm(object, ...)}: Fits HMM to BeadStudioSetList
      object.  Additional arguments can be passed to
      \code{hmmBeadStudioSetList}.}

    \item{}{\code{length(x)}:

      Returns the number of elements in the list object.

    }

  }
}

\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{BeadStudioSetList}}
}
\examples{
new("BeadStudioSetList")
}
\keyword{methods}
