\name{HmmOptions-class}
\alias{[,HmmOptions-method}
\alias{annotation,HmmOptions-method}
\alias{annotation<-,HmmOptions,ANY-method}
\alias{db,HmmOptions-method}
\alias{hmm,HmmOptions,HmmParameter-method}
\alias{initialize,HmmOptions-method}
\alias{notes,HmmOptions-method}
\alias{show,HmmOptions-method}
\alias{calculateDistance,HmmOptions-method}
\alias{physicalDistance,HmmOptions-method}
\alias{calls.ICE,HmmOptions-method}
\alias{calls.ICE<-,HmmOptions,logical-method}
\alias{copyNumber.ICE,HmmOptions-method}
\alias{copyNumber.ICE<-,HmmOptions,logical-method}
\alias{copyNumber.location,HmmOptions-method}
\alias{copyNumber.location,HmmOptions,numeric-method}
\alias{probHomCall,HmmOptions-method}
\alias{probHomCall<-,HmmOptions,numeric-method}
\alias{snpset,HmmOptions-method}
\alias{snpset<-,HmmOptions,SnpLevelSet-method}
\alias{states,HmmOptions-method}
\alias{states<-,HmmOptions,character-method}
\alias{term5,HmmOptions-method}
\alias{term5<-,HmmOptions,numeric-method}
\docType{class}
\alias{HmmOptions-class}
\title{Class "HmmOptions"}
\description{Options for fitting the hidden Markov model.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HmmOptions", snpset, states, copyNumber.location, copyNumber.scale, copyNumber.ICE, calls.ICE, probHomCall, term5)}.
}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{Object of class \code{"character"}.  Hidden
      states of the HMM.}
    \item{\code{snpset}:}{\code{"SnpLevelSet"} instance. }
    \item{\code{copyNumber.location}:}{Object of class
      \code{"numeric"}.  Ignored when snpset is an instance of
      \code{SnpCallSet}. See \code{\link{copyNumber.location}}}
    
    \item{\code{copyNumber.scale}:}{Object of class
      \code{"NULLOrNumeric"}.  This slot is not currently used.  If one
      has standard errors for the copy number estimates, the inverse of
      the standard errors should be stored in the \code{cnConfidence}
      \code{assayData} element of the snpset object.}
    
    \item{\code{copyNumber.ICE}:}{Object of class \code{"logical"}. See \code{\link{copyNumber.ICE}}.  }
    
    \item{\code{calls.ICE}:}{Object of class \code{"logical"}.  See \code{\link{calls.ICE}}}
    \item{\code{probHomCall}:}{Object of class \code{"numeric"}. The
    probability of a homozygous genotype call corresponding to each of
    the states. When \code{calls.ICE} is TRUE, only the first two
    elements will be used.}
    \item{\code{term5}:}{Object of class
    \code{"NULLOrNumeric"}. User-specified probabilities of f(HOM |
    genotype call, state).  Ignored when
    \code{calls.ICE} is \code{FALSE}.  For more details, see Equation
    2.7 in the manuscript that this package references.}
  }
}
\section{Methods}{
  \describe{

    \item{[}{\code{signature(x = "HmmOptions")}: }
    \item{annotation}{\code{signature(object = "HmmOptions")}: }
    \item{annotation<-}{\code{signature(object = "HmmOptions", value = "ANY")}: }
    \item{calculateDistance}{\code{signature(object = "HmmOptions")}: }
    \item{calls.ICE}{\code{signature(object = "HmmOptions")}: }
    \item{calls.ICE<-}{\code{signature(object = "HmmOptions", value="logical")}: }    
    \item{copyNumber.ICE}{\code{signature(object = "HmmOptions")}: }
    \item{copyNumber.ICE<-}{\code{signature(object = "HmmOptions", value="logical")}: }        
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
