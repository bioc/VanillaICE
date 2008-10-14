\name{HmmParameter-class}
\docType{class}
\alias{HmmParameter-class}
\alias{[,HmmParameter-method}
\alias{emission,HmmParameter-method}
\alias{emission<-,HmmParameter,array-method}
\alias{genomicDistance,HmmParameter-method}
\alias{genomicDistance<-,HmmParameter,numeric-method}
\alias{initialize,HmmParameter-method}
\alias{Pi,HmmParameter-method}
\alias{show,HmmParameter-method}
\alias{states,HmmParameter-method}
\alias{states<-,HmmParameter,character-method}
\alias{tau,HmmParameter-method}
\alias{transitionScale,HmmParameter-method}
\alias{transitionScale<-,HmmParameter,matrix-method}

\title{Class "HmmParameter" }
\description{Parameters needed for the Viterbi algorithm, including
  emission probabilities, the genomic distance, and the initial state probabilities}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HmmParameter", states, initialStateProbability, emission, genomicDistance, transitionScale)}.
}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{Object of class \code{"character"} }
    \item{\code{initialStateProbability}:}{Object of class \code{"numeric"} }
    \item{\code{emission}:}{Object of class \code{"array"} }
    \item{\code{genomicDistance}:}{Object of class \code{"numeric"} }
    \item{\code{transitionScale}:}{Object of class \code{"matrix"} }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "HmmParameter")}:  }
    \item{emission}{\code{signature(object = "HmmParameter")}:  }
    \item{emission<-}{\code{signature(object = "HmmParameter", value = "array")}:  }
    \item{genomicDistance}{\code{signature(object = "HmmParameter")}:  }
    \item{genomicDistance<-}{\code{signature(object = "HmmParameter", value = "numeric")}:  }
    \item{hmm}{\code{signature(object = "HmmOptions", params = "HmmParameter")}:  }
    \item{initialize}{\code{signature(.Object = "HmmParameter")}:  }
    \item{pi}{\code{signature(object = "HmmParameter")}:  }
    \item{show}{\code{signature(object = "HmmParameter")}:  }
    \item{states}{\code{signature(object = "HmmParameter")}:  }
    \item{states<-}{\code{signature(object = "HmmParameter", value = "character")}:  }
    \item{tau}{\code{signature(object = "HmmParameter")}:  }
    \item{transitionScale}{\code{signature(object = "HmmParameter")}:  }
    \item{transitionScale<-}{\code{signature(object = "HmmParameter", value = "matrix")}:  }
	 }
}
\author{RS}
\note{See HmmParameter vignette}
\seealso{
  \code{\link{HmmOptions-class}}
}
\examples{
showClass("HmmParameter")
}
\keyword{classes}
