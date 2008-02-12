\name{ParHSet-class}
\docType{class}
\alias{ParHSet-class}
\alias{getPar,ParHSet,SnpLevelSet-method}
\alias{initialize,ParHSet-method}
\alias{plotSnp,ParHSet,SnpLevelSet-method}
\alias{show,ParHSet-method}
\title{Graphical parameters for plotting}
\description{Generating an instance of this class creates a set of
  default graphical parameters useful for plotting the SNP-level data
  (genotype and copy number) as well as the predicted states from the
  hidden Markov model}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParHSet", border.prediction, legend.prediction, legend.location.prediction, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}:  a list of
      graphical parameters. See \code{\link[graphics]{par}}}
  }
}
\section{Extends}{
Class \code{\link[SNPchip]{ParESet-class}}, directly.
}
\section{Methods}{
  \describe{
    \item{getPar}{\code{signature(.Object = "ParHSet")} }
    \item{initialize}{\code{signature(.Object = "ParHSet")}}
    \item{plotSnp}{\code{signature(object = "ParHSet", snpset =
	"SnpLevelSet")}: See also \code{\link[SNPchip]{plotSnp}}}
    \item{show}{\code{signature(object = "ParHSet")}}
	 }
}
\author{R. Scharpf}
\seealso{
  \code{\link[SNPchip]{plotSnp}}
}
\examples{
showClass("ParHSet")
}
\keyword{classes}
