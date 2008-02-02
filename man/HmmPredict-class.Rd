\name{HmmPredict-class}
\docType{class}
\alias{HmmPredict-class}
\alias{SnpClass,HmmPredict-method}
\alias{[,HmmPredict-method}
\alias{breakpoints}
\alias{breakpoints,HmmPredict-method}
\alias{featureNames,HmmPredict-method}
\alias{getPar,HmmPredict-method}
\alias{initialize,HmmPredict-method}
\alias{predictions}
\alias{predictions,HmmPredict-method}
\alias{sampleNames,HmmPredict-method}
\alias{show,HmmPredict-method}
\alias{states,HmmPredict-method}
\title{Class for the hidden Markov model results}
\description{Contains the predictions for the hidden states and a list
  of breakpoint matrices.  Each element in the list corresponds to the
  breakpoints identified in one sample.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HmmPredict", states, predictions, breakpoints, SnpClass)}.
}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{Object of class \code{"character"}. Character
      string of the hidden states.}
    \item{\code{predictions}:}{Object of class \code{"matrix"}. Matrix
      of predictions.}
    \item{\code{breakpoints}:}{Object of class \code{"list"}.  List of
      matrices.  Each element is the matrix of breakpoints for one sample.}
    \item{\code{SnpClass}:}{Object of class \code{"character"}. The
      class of SNP data (e.g., \code{oligoSnpSet}) }
  }
}
\section{Methods}{
  \describe{
    \item{SnpClass}{\code{signature(object = "HmmPredict")}: Accessor
      for class of SNP data}
    \item{[}{\code{signature(x = "HmmPredict")}: Subsetting the predictions}
    \item{breakpoints}{\code{signature(object = "HmmPredict")}: Accessor
    for the list of breakpoints}
    \item{featureNames}{\code{signature(object = "HmmPredict")}:
      Character string of SNP identifiers}
    \item{getPar}{\code{signature(object = "HmmPredict")}: Graphical parameters.}
    \item{initialize}{\code{signature(.Object = "HmmPredict")}  }
    \item{predictions}{\code{signature(object = "HmmPredict")}: Accessor
    for the predictions}
    \item{sampleNames}{\code{signature(object = "HmmPredict")}: Sample names}
    \item{show}{\code{signature(object = "HmmPredict")}}
    \item{states}{\code{signature(object = "HmmPredict")}: Accessor for
      hidden states}
	 }
}
\author{R. Scharpf}

\seealso{
  \code{\link{hmm}}
}
\examples{
showClass("HmmPredict")
}
\keyword{classes}
