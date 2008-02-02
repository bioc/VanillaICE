\name{ParHmmSnpSet-class}
\docType{class}
\alias{ParHmmSnpSet-class}
\alias{initialize,ParHmmSnpSet-method} 
\alias{plotSnp,ParHmmSnpSet,SnpLevelSet-method}
\title{Graphical parameters for plotting}
\description{Graphical parameters for plotting when the data object is
  an object of class \code{oligoSnpSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParHmmSnpSet", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}. A list of
      graphical parameters }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ParHSet}"}, directly.
Class \code{"\linkS4class{ParESet}"}, by class "ParHSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParHmmSnpSet")}}
    \item{plotSnp}{\code{signature(object = "ParHmmSnpSet", snpset =
	"SnpLevelSet")} See also \code{\link[SNPchip]{plotSnp}}}
  }
}
\author{R. Scharpf}

\seealso{
  \code{\link[SNPchip]{plotSnp}}, 
  \code{\link[oligoClasses]{oligoSnpSet-class}}, \code{\link[graphics]{par}}
}
\examples{
showClass("ParHmmSnpSet")
}
\keyword{classes}
