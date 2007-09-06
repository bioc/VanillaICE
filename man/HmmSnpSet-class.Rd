\name{HmmSnpSet-class}
\docType{class}
\alias{HmmSnpSet-class}
\alias{[,HmmSnpSet-method}
\alias{[<-,HmmSnpSet-method}
\alias{bind,HmmSnpSet,HmmSnpSet-method}
\alias{bind,NULL,HmmSnpSet-method}
\alias{bind,HmmSnpSet,NULL-method}
%\alias{breakpoints,list-method}
\alias{breakpoints,HmmSnpSet-method}
\alias{calculateEmissionProbability,HmmSnpSet-method}
\alias{coerce,oligoSnpSet,HmmSnpSet-method}
\alias{getCopyNumberLocation,HmmSnpSet-method}
\alias{getCopyNumberScale,HmmSnpSet-method}
\alias{getInitialStateProbability,HmmSnpSet-method}
\alias{getLohProbability,HmmSnpSet-method}
\alias{getTransitionProbability,HmmSnpSet-method}
\alias{hmm,HmmSnpSet-method}
\alias{ice2VanillaICE}
\alias{ice2VanillaICE,HmmSnpSet-method}
\alias{initialize,HmmSnpSet-method}
\alias{locationCopyNumber,HmmSnpSet-method}
\alias{locationCopyNumber<-,HmmSnpSet-method}
\alias{pCHOM,HmmSnpSet-method}
\alias{pCHOM<-,HmmSnpSet-method}
\alias{predictions,HmmSnpSet-method}
\alias{predictions<-,HmmSnpSet,matrix-method}
\alias{scaleCopyNumber,HmmSnpSet-method}
\alias{scaleCopyNumber<-,HmmSnpSet-method}
\alias{show,HmmSnpSet-method}
\alias{.unsplitHmm}
\alias{.unsplitHmm,list,AnnotatedDataFrame-method}

\title{HmmSnpSet}
\description{}
\section{Objects from the Class}{
  
Objects can be created by calls of the form \code{new("HmmSnpSet",
assayData, featureData, phenoData, experimentData, annotation, calls,
callsConfidence, copyNumber, cnConfidence, predictions, scaleCopyNumber,
probabilityCalledHomIfNormal, probabilityCalledHetIfLoh, log2,
locationCopyNumber, stateNames, ...)}.

}

\section{Slots}{
	 \describe{
    \item{\code{locationCopyNumber}:}{Object of class \code{"numeric"}  }
    \item{\code{scaleCopyNumber}:}{Object of class \code{"NullOrMatrix"}  }
    \item{\code{log2}:}{Object of class \code{"logical"}  }
    \item{\code{predictions}:}{Object of class \code{"matrix"}  }
    \item{\code{distance}:}{Object of class \code{"logical"}  }
    \item{\code{initialStateProbability}:}{Object of class \code{"NullOrNumeric"}  }
    \item{\code{transitionProbability}:}{Object of class \code{"NullOrNumericOrMatrix"}  }
    \item{\code{emissionProbability}:}{Object of class \code{"NullOrMatrix"}  }
    \item{\code{likelihood}:}{Object of class \code{"NullOrMatrix"}  }
    \item{\code{chromosomeArm}:}{Object of class \code{"character"}  }
    \item{\code{callsIce}:}{Object of class \code{"logical"}  }    
    \item{\code{copyNumberIce}:}{Object of class \code{"logical"}  }    
    \item{\code{stateNames}:}{Object of class \code{"character"}  }
    \item{\code{chromosomeAnnotation}:}{Object of class \code{"data.frame"}  }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
    \item{\code{probabilityCalledHomIfNormal}:}{Object of class \code{"numeric"}  }
    \item{\code{probabilityCalledHetIfLoh}:}{Object of class \code{"numeric"}  }
  }
}
\section{Extends}{
Class \code{\link{HmmSnpCopyNumberSet-class}}, directly.
Class \code{\link{HmmSnpCallSet-class}}, directly.
Class \code{\link{hSet-class}}, by class "HmmSnpCopyNumberSet", distance 2.
Class \code{\link[SNPchip]{AnnotatedSnpSet-class}}, by class "HmmSnpCopyNumberSet", distance 3.
Class \code{\link[SNPchip]{AnnotatedSnpCallSet-class}}, by class "HmmSnpCopyNumberSet", distance 4.
Class \code{\link[SNPchip]{AnnotatedSnpCopyNumberSet-class}}, by class "HmmSnpCopyNumberSet", distance 4, with explicit coerce.
Class \code{\link[oligo]{SnpCallSet-class}}, by class "HmmSnpCopyNumberSet", distance 5.
Class \code{\link[Biobase]{class.eSet}}, by class "HmmSnpCopyNumberSet", distance 6.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "HmmSnpCopyNumberSet", distance 7.
Class \code{\link[Biobase]{Versioned-class}}, by class "HmmSnpCopyNumberSet", distance 8.
Class \code{\link[oligo]{SnpCopyNumberSet-class}}, by class "HmmSnpCopyNumberSet", distance 5, with explicit coerce.
Class \code{\link[Biobase]{class.eSet}}, by class "HmmSnpCopyNumberSet", distance 6, with explicit coerce.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "HmmSnpCopyNumberSet", distance 7, with explicit coerce.
Class \code{\link[Biobase]{Versioned-class}}, by class "HmmSnpCopyNumberSet", distance 8, with explicit coerce.
Class \code{\link{hSet-class}}, by class "HmmSnpCallSet", distance 2.
Class \code{\link[SNPchip]{AnnotatedSnpSet-class}}, by class "HmmSnpCallSet", distance 3.
Class \code{\link[SNPchip]{AnnotatedSnpCallSet-class}}, by class "HmmSnpCallSet", distance 4.
Class \code{\link[SNPchip]{AnnotatedSnpCopyNumberSet-class}}, by class "HmmSnpCallSet", distance 4, with explicit coerce.
Class \code{\link[oligo]{SnpCallSet-class}}, by class "HmmSnpCallSet", distance 5.
Class \code{\link[Biobase]{class.eSet}}, by class "HmmSnpCallSet", distance 6.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "HmmSnpCallSet", distance 7.
Class \code{\link[Biobase]{Versioned-class}}, by class "HmmSnpCallSet", distance 8.
Class \code{\link[oligo]{SnpCopyNumberSet-class}}, by class "HmmSnpCallSet", distance 5, with explicit coerce.
Class \code{\link[Biobase]{class.eSet}}, by class "HmmSnpCallSet", distance 6, with explicit coerce.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "HmmSnpCallSet", distance 7, with explicit coerce.
Class \code{\link[Biobase]{Versioned-class}}, by class "HmmSnpCallSet", distance 8, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "HmmSnpSet")}: ... }
    \item{[<-}{\code{signature(x = "HmmSnpSet")}: ... }
    \item{breakpoints}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{calculateEmissionProbability}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{coerce}{\code{signature(from = "AnnotatedSnpSet", to =
	"HmmSnpSet")}: ... }
    \item{callsIce}{\code{signature(object = "HmmSnpSet")}:
      ... }
    \item{copyNumberIce}{\code{signature(object = "HmmSnpSet")}: ... }    
    \item{getCopyNumberLocation}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{getLohProbability}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{getTransitionProbability}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{hmm}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{ice2VanillaICE}{\code{signature(.Object = "HmmSnpSet")}:
	 Resets the package attribute from ICE (DEPRECATED package name)
	 to VanillaICE}
    \item{initialize}{\code{signature(.Object = "HmmSnpSet")}: ... }
    \item{locationCopyNumber}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{locationCopyNumber<-}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{plotHmm}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{predictions}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{predictions<-}{\code{signature(object = "HmmSnpSet", value = "matrix")}: ... }
    \item{scaleCopyNumber}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{scaleCopyNumber<-}{\code{signature(object = "HmmSnpSet")}: ... }
    \item{show}{\code{signature(object = "HmmSnpSet")}: ... }
	 }
}
\references{}
\author{Robert Scharpf}
\note{ }

\seealso{
  
  \code{\link{HmmSnpCopyNumberSet-class}},
  \code{\link{HmmSnpCallSet-class}}
  \code{\link{hSet-class}}
}

\examples{
showClass("HmmSnpSet")
}
\keyword{classes}
