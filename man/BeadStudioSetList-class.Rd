\name{BeadStudioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{BeadStudioSetList-class}
\alias{[,BeadStudioSetList-method}
\alias{[[,BeadStudioSetList-method}
\alias{$,BeadStudioSetList-method}
\alias{$<-,BeadStudioSetList-method}
\alias{hmm,BeadStudioSetList-method}
\alias{length,BeadStudioSetList-method}

\title{BeadStudioSetList class}

\description{
  Container for log R ratios and B allele frequencies stored by chromosome.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BeadStudioSetList", assayDataList, logRRatio, BAF, featureDataList, chromosome, phenoData, annotation, genomeBuild, ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{assayDataList}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureDataList}:}{Object of class \code{"list"} ~~ }
    \item{\code{chromosome}:}{Object of class \code{"integer"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{genomeBuild}:}{Object of class \code{"character"} ~~ }
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "BeadStudioSetList")}: ... }
    \item{}{\code{object[[i]]}: Returns an object of class
      \code{BeadStudioSet} from the ith element in \code{object}.}
    \item{$}{\code{signature(x = "BeadStudioSetList")}: ... }
    \item{}{\code{object$name <- value}: assign vector \code{value} to
      variable name \code{name} of \code{phenoData}}
    \item{}{\code{hmm(object, ...)}: Fits HMM to BeadStudioSetList
      object.  Additional arguments can be passed to
      \code{hmmBeadStudioSetList}.}
    \item{length}{\code{signature(x = "BeadStudioSetList")}: ... }
  }
}

\author{
R. Scharpf
}

\seealso{
  \code{\link{BeadStudioSetList}}
}
\examples{
new("BeadStudioSetList")
}
\keyword{classes}
