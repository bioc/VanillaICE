\name{ParHmmSnpSet-class}
\docType{class}
\alias{ParHmmSnpSet-class}
\alias{initialize,ParHmmSnpSet-method}
\alias{plotSnp,ParHmmSnpSet,HmmSnpSet-method}
\title{Graphical parameters }
\description{Graphical parameters for plotting object of class \code{HmmSnpSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParHmmSnpSet", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}}
  }
}
\section{Extends}{
Class \code{\link{ParHSet-class}}, directly.
Class \code{\link[SNPchip]{ParESet-class}}, by class "ParHSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParHmmSnpSet")}: ... }
  }
}
\references{}
\author{}
\note{}

\seealso{
  \code{\link[SNPchip]{plotSnp}}, \code{\link[SNPchip]{ParESet-class}},
  \code{\link{ParHSet-class}}
}
\examples{
showClass("ParHmmSnpSet")
}
\keyword{classes}
