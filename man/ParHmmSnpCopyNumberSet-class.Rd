\name{ParHmmSnpCopyNumberSet-class}
\docType{class}
\alias{ParHmmSnpCopyNumberSet-class}
\alias{initialize,ParHmmSnpCopyNumberSet-method}
\alias{plotSnp,ParHmmSnpCopyNumberSet,SnpLevelSet-method}
\title{Graphical parameters for plotting}
\description{Graphical parameters for plotting data of class
  \code{SnpCopyNumberSet} and the corresponding hidden states predicted
  from the hidden Markov model}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParHmmSnpCopyNumberSet", border.prediction, legend.prediction, legend.location.prediction, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}. A list of
      graphical parameters. See \code{par}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ParHSet}"}, directly.
Class \code{"\linkS4class{ParESet}"}, by class "ParHSet", distance 2.
}
\section{Methods}{
  \item{initialize}{\code{signature(.Object = "ParHmmSnpSet")}}
  \item{plotSnp}{\code{signature(object = "ParHmmSnpCopyNumberSet", snpset =
	"SnpLevelSet")} See also \code{\link[SNPchip]{plotSnp}}}
}
\author{R. Scharpf}
\seealso{
  \code{\link[SNPchip]{plotSnp}}, 
  \code{\link[oligoClasses]{SnpCopyNumberSet-class}}, \code{\link[graphics]{par}}
}
\examples{
showClass("ParHmmSnpCopyNumberSet")
}
\keyword{classes}
