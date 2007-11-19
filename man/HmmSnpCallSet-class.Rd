\name{HmmSnpCallSet-class}
\docType{class}
\alias{HmmSnpCallSet-class}
\alias{[,HmmSnpCallSet-method}
\alias{bind,HmmSnpCallSet,HmmSnpCallSet-method}
\alias{bind,NULL,HmmSnpCallSet-method}
\alias{bind,HmmSnpCallSet,NULL-method}
\alias{calculateEmissionProbability,HmmSnpCallSet-method}
\alias{callsIce,HmmSnpCallSet-method}
\alias{callsIce<-,HmmSnpCallSet,logical-method}
\alias{coerce,oligoSnpSet,HmmSnpCallSet-method}
\alias{coerce,HmmSnpSet,HmmSnpCallSet-method}
\alias{coerce,SnpCallSet,HmmSnpCallSet-method}
\alias{getInitialStateProbability,HmmSnpCallSet-method}
\alias{getLohProbability,HmmSnpCallSet-method}
\alias{getTransitionProbability,HmmSnpCallSet-method}
\alias{hmm,HmmSnpCallSet-method}
\alias{initialize,HmmSnpCallSet-method}
\alias{isMissing,HmmSnpCallSet-method}
\alias{pCHOM,HmmSnpCallSet-method}
\alias{pCHOM<-,HmmSnpCallSet-method}

\alias{show,HmmSnpCallSet-method}

\title{Class HmmSnpCallSet}

\description{A class of objects for containing summaries of SNP-level,
  as well as predictions from the HMM model.  In particular, the
  HmmSnpCallSet class is an extension of SnpCallSet, requiring only
  SNP-level summaries of genotype calls. }

\section{Objects from the Class}{
  
Objects can be created by calls of the form \code{new("HmmSnpCallSet",
assayData, phenoData, featureData, experimentData, annotation,
chromosomeAnnotation, calls, callsConfidence, copyNumber, cnConfidence,
...)}.

}
\section{Slots}{
	 \describe{
    \item{\code{probabilityCalledHomIfNormal}:}{Object of class \code{"numeric"}  }
    \item{\code{probabilityCalledHetIfLoh}:}{Object of class \code{"numeric"}  }
    \item{\code{predictions}:}{Object of class \code{"matrix"}  }
    \item{\code{distance}:}{Object of class \code{"logical"}  }
    \item{\code{initialStateProbability}:}{Object of class \code{"NullOrNumeric"}  }
    \item{\code{transitionProbability}:}{Object of class \code{"NullOrNumericOrMatrix"}  }
    \item{\code{emissionProbability}:}{Object of class \code{"NullOrMatrix"}  }
    \item{\code{likelihood}:}{Object of class \code{"NullOrMatrix"}  }
    \item{\code{chromosomeArm}:}{Object of class \code{"character"}  }
    \item{\code{callsIce}:}{Object of class \code{"logical"}}
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
Class \code{\link{hSet-class}}, directly.
Class \code{\link[SNPchip]{oligoSnpSet-class}}, by class "hSet", distance 2.
Class \code{\link[SNPchip]{SnpCallSet-class}}, by class "hSet", distance 3.
Class \code{\link[SNPchip]{SnpCopyNumberSet-class}}, by class "hSet", distance 3, with explicit coerce.
Class \code{\link[Biobase]{SnpCallSet-class}}, by class "hSet", distance 4.
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
    \item{[}{\code{signature(x = "HmmSnpCallSet")} }
    \item{calculateEmissionProbability}{\code{signature(object =
	"HmmSnpCallSet")}}
    \item{callsIce}{\code{signature(object = "HmmSnpCallSet")} }    
    \item{coerce}{\code{signature(from = "AnnotatedSnpSet", to = "HmmSnpCallSet")} }
    \item{getLohProbability}{\code{signature(object = "HmmSnpCallSet")} }
    \item{getTransitionProbability}{\code{signature(object = "HmmSnpCallSet")} }
    \item{hmm}{\code{signature(object = "HmmSnpCallSet")} }
    \item{isMissing}{\code{signature(object = "HmmSnpCallSet")} }
    \item{plotHmm}{\code{signature(object = "HmmSnpCallSet")} }
    \item{predictions}{\code{signature(object = "HmmSnpCallSet")} }
    \item{predictions<-}{\code{signature(object = "HmmSnpCallSet", value = "ANY")} }
    \item{show}{\code{signature(object = "HmmSnpCallSet")} }
	 }
}
\references{}
\author{Robert Scharpf}
\note{}

\seealso{

  \code{\link{HmmSnpSet-class}},
  \code{\link{HmmSnpCopyNumberSet-class}},
  \code{\link{hSet-class}}
  
}

\examples{
showClass("HmmSnpCallSet")
}
\keyword{classes}
