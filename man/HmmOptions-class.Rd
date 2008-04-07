\name{HmmOptions-class}
\docType{class}
\alias{HmmOptions-class}
\alias{[,HmmOptions-method}
\alias{annotation,HmmOptions-method}
\alias{annotation<-,HmmOptions,ANY-method}
\alias{calculateDistance,HmmOptions-method}
\alias{db,HmmOptions-method}
\alias{hmm,HmmOptions,HmmParameter-method}
\alias{initialize,HmmOptions-method}
\alias{notes,HmmOptions-method}
\alias{show,HmmOptions-method}
\alias{states,HmmOptions-method}
\alias{states<-,HmmOptions,character-method}

\title{Class "HmmOptions"}
\description{Options for fitting the hidden Markov model.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HmmOptions", snpset, states, copyNumber.location, copyNumber.scale, copyNumber.ICE, calls.ICE, probHomCall, probWrongCall)}.
}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{Object of class \code{"character"} }
    \item{\code{snpset}:}{Object of class \code{"SnpLevelSet"} }
    \item{\code{copyNumber.location}:}{Object of class \code{"numeric"}}
    \item{\code{copyNumber.scale}:}{Object of class \code{"NULLOrNumeric"}}
    \item{\code{copyNumber.ICE}:}{Object of class \code{"logical"}}
    \item{\code{calls.ICE}:}{Object of class \code{"logical"}}
    \item{\code{probHomCall}:}{Object of class \code{"numeric"}}
    \item{\code{probWrongCall}:}{Object of class \code{"NULLOrNumeric"}}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "HmmOptions")}: }
    \item{annotation}{\code{signature(object = "HmmOptions")}: }
    \item{annotation<-}{\code{signature(object = "HmmOptions", value = "ANY")}: }
    \item{calculateDistance}{\code{signature(object = "HmmOptions")}: }
    \item{db}{\code{signature(object = "HmmOptions")}: }
    \item{hmm}{\code{signature(object = "HmmOptions", params = "HmmParameter")}:}
    \item{initialize}{\code{signature(.Object = "HmmOptions")}: }
    \item{notes}{\code{signature(object = "HmmOptions")}: }
    \item{show}{\code{signature(object = "HmmOptions")}:  }
    \item{states}{\code{signature(object = "HmmOptions")}:  }
    \item{states<-}{\code{signature(object = "HmmOptions", value = "character")}: }
	 }
}
\author{R. Scharpf}
\seealso{
  \code{\link{HmmParameter-class}}
}
\examples{
showClass("HmmOptions")
}
\keyword{classes}
