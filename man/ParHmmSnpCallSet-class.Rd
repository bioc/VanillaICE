\name{ParHmmSnpCallSet-class}
\docType{class}
\alias{ParHmmSnpCallSet-class}
\alias{initialize,ParHmmSnpCallSet-method} 
\title{Graphical parameters for plotting}
\description{Graphical parameters for plotting \code{SnpCallSet} and the
corresponding predicted states from the hidden Markov model}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParHmmSnpCallSet", border.prediction, legend.prediction, legend.location.prediction, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}. List of
      graphical parameters.  See \code{par}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ParHSet}"}, directly.
Class \code{"\linkS4class{ParESet}"}, by class "ParHSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParHmmSnpSet")}}
    }
}
\author{R. Scharpf}
\seealso{
  \code{\link[SNPchip]{plotSnp}}
  \code{\link[oligoClasses]{SnpCallSet-class}}, \code{\link[graphics]{par}}
}
\examples{
showClass("ParHmmSnpCallSet")
}
\keyword{classes}
