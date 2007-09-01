\name{ParHmmSnpCallSet-class}
\docType{class}
\alias{ParHmmSnpCallSet-class}
\alias{initialize,ParHmmSnpCallSet-method}

\title{Class ParHmmSnpCallSet}
\description{Graphical parameters for objects of class \code{HmmSnpCallSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParHmmSnpCallSet", ...)}.
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
    \item{initialize}{\code{signature(.Object = "ParHmmSnpCallSet")}}
  }
}
\references{}
\author{}
\note{}

\seealso{
  \code{\link{HmmSnpCallSet-class}}
}
\examples{
showClass("ParHmmSnpCallSet")
snpPar(new("ParHmmSnpCallSet"))
}
\keyword{classes}
