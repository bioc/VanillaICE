\name{HmmSnpCopyNumberSet-class}
\docType{class}
\alias{HmmSnpCopyNumberSet-class}
\alias{[,HmmSnpCopyNumberSet-method}
\alias{bind,HmmSnpCopyNumberSet,HmmSnpCopyNumberSet-method}
\alias{bind,NULL,HmmSnpCopyNumberSet-method}
\alias{bind,HmmSnpCopyNumberSet,NULL-method}
\alias{calculateEmissionProbability,HmmSnpCopyNumberSet-method}
\alias{coerce,oligoSnpSet,HmmSnpCopyNumberSet-method}
\alias{coerce,HmmSnpSet,HmmSnpCopyNumberSet-method}
\alias{coerce,SnpCopyNumberSet,HmmSnpCopyNumberSet-method}
\alias{copyNumberIce,HmmSnpCopyNumberSet-method}
\alias{copyNumberIce<-,HmmSnpCopyNumberSet,logical-method}
\alias{getCopyNumberLocation,HmmSnpCopyNumberSet-method}
\alias{getCopyNumberScale,HmmSnpCopyNumberSet-method}
\alias{getInitialStateProbability,HmmSnpCopyNumberSet-method}
\alias{getTransitionProbability,HmmSnpCopyNumberSet-method}
\alias{hmm,HmmSnpCopyNumberSet-method}
\alias{initialize,HmmSnpCopyNumberSet-method}
\alias{locationCopyNumber,HmmSnpCopyNumberSet-method}
\alias{locationCopyNumber<-,HmmSnpCopyNumberSet-method}
\alias{predictions,HmmSnpCopyNumberSet-method}
\alias{predictions<-,HmmSnpCopyNumberSet,ANY-method}
\alias{scaleCopyNumber,HmmSnpCopyNumberSet-method}
\alias{scaleCopyNumber<-,HmmSnpCopyNumberSet-method}
\alias{show,HmmSnpCopyNumberSet-method}

\title{HmmSnpCopyNumberSet}
\description{}
\section{Objects from the Class}{
  
Objects can be created by calls of the form
\code{new("HmmSnpCopyNumberSet", assayData, phenoData, featureData,
experimentData, annotation, chromosomeAnnotation, calls,
callsConfidence, copyNumber, cnConfidence, ...)}.


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
    \item{\code{stateNames}:}{Object of class \code{"character"}  }
    \item{\code{chromosomeAnnotation}:}{Object of class \code{"data.frame"}  }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
    \item{\code{copyNumberIce}:}{Object of class \code{"logical"} }          
  }
}
\section{Extends}{
Class \code{\link{hSet-class}}, directly.
Class \code{\link[SNPchip]{AnnotatedSnpSet-class}}, by class "hSet", distance 2.
Class \code{\link[SNPchip]{AnnotatedSnpCallSet-class}}, by class "hSet", distance 3.
Class \code{\link[SNPchip]{AnnotatedSnpCopyNumberSet-class}}, by class "hSet", distance 3, with explicit coerce.
Class \code{\link[oligo]{SnpCallSet-class}}, by class "hSet", distance 4.
Class \code{\link[Biobase]{class.eSet}}, by class "hSet", distance 5.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "hSet", distance 6.
Class \code{\link[Biobase]{Versioned-class}}, by class "hSet", distance 7.
Class \code{\link[oligo]{SnpCopyNumberSet-class}}, by class "hSet", distance 4, with explicit coerce.
Class \code{\link[Biobase]{class.eSet}}, by class "hSet", distance 5, with explicit coerce.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "hSet", distance 6, with explicit coerce.
Class \code{\link[Biobase]{Versioned-class}}, by class "hSet", distance 7, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "HmmSnpCopyNumberSet")} }
    \item{calculateEmissionProbability}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{coerce}{\code{signature(from = "AnnotatedSnpSet", to =
	"HmmSnpCopyNumberSet")}}
    \item{copyNumberIce}{\code{signature(object = "HmmSnpCopyNumberSet")} }    
    \item{getCopyNumberLocation}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{getTransitionProbability}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{hmm}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{locationCopyNumber}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{locationCopyNumber<-}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{plotHmm}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{predictions}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{predictions<-}{\code{signature(object = "HmmSnpCopyNumberSet", value = "ANY")} }
    \item{scaleCopyNumber}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{scaleCopyNumber<-}{\code{signature(object = "HmmSnpCopyNumberSet")} }
    \item{show}{\code{signature(object = "HmmSnpCopyNumberSet")} }
	 }
}
\references{}
\author{Robert Scharpf}
\note{}

\seealso{

    \code{\link{HmmSnpSet-class}},
    \code{\link{HmmSnpCallSet-class}}
    \code{\link{hSet-class}}    
    
}
\examples{
showClass("HmmSnpCopyNumberSet")
}
\keyword{classes}
