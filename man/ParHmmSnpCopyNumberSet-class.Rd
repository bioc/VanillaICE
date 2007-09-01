\name{ParHmmSnpCopyNumberSet-class}
\docType{class}
\alias{ParHmmSnpCopyNumberSet-class}
\alias{initialize,ParHmmSnpCopyNumberSet-method}
\title{Class ParHmmSnpCopyNumberSet}
\description{Graphical parameters for objects of class \code{HmmSnpCopyNumberSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParHmmSnpCopyNumberSet", ...)}.
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
    \item{initialize}{\code{signature(.Object = "ParHmmSnpCopyNumberSet")}}
  }
}
\references{}
\author{}
\note{}
\seealso{

}
\examples{
showClass("ParHmmSnpCopyNumberSet")
snpPar(new("ParHmmSnpCopyNumberSet"))
}
\keyword{classes}
