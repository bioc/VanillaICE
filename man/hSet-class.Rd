\name{hSet-class}
\docType{class}
\alias{hSet-class}
\alias{breakpoints,hSet-method}
\alias{byArm,hSet-method}
\alias{callsIce,hSet-method}
\alias{callsIce<-,hSet,logical-method}
\alias{copyNumberIce,hSet-method}
\alias{copyNumberIce<-,hSet,logical-method}
\alias{chromosomeArm,hSet-method}
\alias{chromosomeArm<-,hSet-method}
\alias{distance,hSet-method}
\alias{distance<-,hSet-method}
\alias{emissionProbability,hSet-method}
\alias{emissionProbability<-,hSet-method}
\alias{getBreaks}
\alias{getBreaks,hSet-method}
\alias{getPar,hSet-method}
\alias{hmmBySample,hSet-method}
\alias{initialStateProbability,hSet-method}
\alias{initialStateProbability<-,hSet-method}
\alias{isMissing,hSet-method}
\alias{likelihood,hSet-method}
\alias{likelihood<-,hSet-method}
\alias{numberStates,hSet-method}
\alias{plotSnp,hSet-method}
\alias{predictions,hSet-method}
\alias{predictions<-,hSet,ANY-method}
\alias{stateNames,hSet-method}
\alias{stateNames<-,hSet-method}
\alias{transitionProbability,hSet-method}
\alias{transitionProbability<-,hSet-method}
\alias{viterbi,hSet-method}

\title{Virtual class for fitting HMMs to SNP chip data}

\description{A virtual class for containing summaries of SNP-level, as
  well as predictions from the HMM model.  In particular, the
  HmmSnpCopyNumberSet, HmmSnpCallSet, and HmmSnpSet classes each extend
  this class. }

\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
	 \describe{
    \item{\code{distance}:}{Object of class \code{"logical"}  }
    \item{\code{initialStateProbability}:}{Object of class \code{"NullOrNumeric"}  }
    \item{\code{transitionProbability}:}{Object of class \code{"NullOrNumericOrMatrix"}  }
    \item{\code{emissionProbability}:}{Object of class \code{"NullOrMatrix"}  }
    \item{\code{likelihood}:}{Object of class \code{"NullOrMatrix"}  }
    \item{\code{chromosomeArm}:}{Object of class \code{"character"}  }
    \item{\code{stateNames}:}{Object of class \code{"character"}  }
    \item{\code{chromosomeAnnotation}:}{Object of class \code{"data.frame"}  }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}
\section{Extends}{
Class \code{\link[SNPchip]{AnnotatedSnpSet-class}}, directly.
Class \code{\link[SNPchip]{AnnotatedSnpCallSet-class}}, by class "AnnotatedSnpSet", distance 2.
Class \code{\link[SNPchip]{AnnotatedSnpCopyNumberSet-class}}, by class "AnnotatedSnpSet", distance 2, with explicit coerce.
Class \code{\link[oligo]{SnpCallSet-clas}}, by class "AnnotatedSnpSet", distance 3.
Class \code{\link[Biobase]{class.eSet}}, by class "AnnotatedSnpSet", distance 4.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "AnnotatedSnpSet", distance 5.
Class \code{\link[Biobase]{Versioned-class}}, by class "AnnotatedSnpSet", distance 6.
Class \code{\link[oligo]{SnpCopyNumberSet-class}}, by class "AnnotatedSnpSet", distance 3, with explicit coerce.
Class \code{\link[Biobase]{class.eSet}}, by class "AnnotatedSnpSet", distance 4, with explicit coerce.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "AnnotatedSnpSet", distance 5, with explicit coerce.
Class \code{\link[Biobase]{Versioned-class}}, by class "AnnotatedSnpSet", distance 6, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{chromosomeArm}{\code{signature(object = "hSet")} }
    \item{chromosomeArm<-}{\code{signature(object = "hSet")} }
    \item{distance}{\code{signature(object = "hSet")} }
    \item{distance<-}{\code{signature(object = "hSet")} }
    \item{emissionProbability}{\code{signature(object = "hSet")} }
    \item{emissionProbability<-}{\code{signature(object = "hSet")} }
    \item{hmmBySample}{\code{signature(object = "hSet")} }
    \item{initialStateProbability}{\code{signature(object = "hSet")} }
    \item{initialStateProbability<-}{\code{signature(object = "hSet")} }
    \item{likelihood}{\code{signature(object = "hSet")} }
    \item{likelihood<-}{\code{signature(object = "hSet")} }
    \item{numberStates}{\code{signature(object = "hSet")} }
    \item{stateNames}{\code{signature(object = "hSet")} }
    \item{stateNames<-}{\code{signature(object = "hSet")} }
    \item{transitionProbability}{\code{signature(object = "hSet")} }
    \item{transitionProbability<-}{\code{signature(object = "hSet")} }
    \item{viterbi}{\code{signature(object = "hSet")} }
	 }
}
\references{Robert Scharpf}
\author{}
\note{}
\seealso{
  
  \code{\link{HmmSnpSet-class}}, \code{\link{HmmSnpCopyNumberSet-class}},
  \code{\link{HmmSnpCallSet-class}}
  
}

\examples{
showClass("hSet")
}
\keyword{classes}
