\name{HmmParameter-class}
\docType{class}
\alias{HmmParameter-class}
\alias{Beta,HmmParameter-method}
\alias{SnpClass}
\alias{SnpClass,HmmParameter-method}
\alias{[,HmmParameter-method}
\alias{arm,HmmParameter-method}
\alias{chromosome,HmmParameter-method}
\alias{cn.ICE,HmmParameter-method}
\alias{cn.location,HmmParameter-method}
\alias{featureData,HmmParameter-method}
\alias{featureNames,HmmParameter-method}
\alias{gt.ICE,HmmParameter-method}
\alias{hmmOptions}
\alias{hmmOptions,HmmParameter-method}
\alias{initialize,HmmParameter-method}
\alias{pData,HmmParameter-method}
\alias{pi}
\alias{pi,HmmParameter-method}
\alias{sampleNames,HmmParameter-method}
\alias{show,HmmParameter-method}
\alias{states}
\alias{states,HmmParameter-method}
\alias{tau}
\alias{tau,HmmParameter-method}
\alias{tau.scale,HmmParameter-method}

\title{Parameters for fitting a hidden Markov model to SNP chip data}
\description{This class contains all the parameters needed to fit the
  hidden Markov model.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HmmParameter", snpset, states, tau, tau.scale, pi, beta, notes, SCALE, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{hmmOptions}:}{Options for fitting the HMM are stored in
      an object of class \code{HmmOptions}}
    \item{\code{tau}:}{Transition probabilities. Object of class
      \code{"matrix"}.  }
    \item{\code{tau.scale}:}{A matrix to scale the transition
      probabilities.  See also \code{\link{tau.scale}}}
    \item{\code{pi}:}{A numerical vector specifying the log initial
      state probabilities}
    \item{\code{beta}:}{A matrix of emission probabilities for the
      observed data.  Each column corresponds to a sample.}
    \item{\code{featureData}:}{Object of class \code{AnnotatedDataFrame}}
    \item{\code{notes}:}{Object of class \code{character}}
  }
}
\section{Methods}{
  \describe{
    \item{Beta}{\code{signature(object = "HmmParameter")} Accessor for
      the emission probabilities.  See also \code{\link{Beta}}}
    \item{SnpClass}{\code{signature(object = "HmmParameter")} A
      character string specifying one of the data classes for high
      throughput SNP data.}
    \item{[}{\code{signature(x = "HmmParameter")} For subsetting.}
    \item{arm}{\code{signature(object = "HmmParameter")} See also \code{\link{arm}}}
    \item{chromosome}{\code{signature(object = "HmmParameter")}:
      Accessor for chromosome information}
    \item{cn.ICE}{\code{signature(object = "HmmParameter")}: See also \code{\link{cn.ICE}}}
    \item{cn.location}{\code{signature(object = "HmmParameter")}: See
      also \code{\link{cn.location}}}
    \item{featureData}{\code{signature(object = "HmmParameter")}:
      Accessor for the \code{\link[Biobase]{featureData}}}
    \item{featureNames}{\code{signature(object = "HmmParameter")}:
      Accessor for the \code{featureNames}}
    \item{gt.ICE}{\code{signature(object = "HmmParameter")}: See also \code{\link{gt.ICE}}}
    \item{hmmOptions}{\code{signature(object = "HmmParameter")}:
      Accessor for the options used to fit the HMM.  See also \code{\link{HmmOptions-class}}}
    \item{initialize}{\code{signature(.Object = "HmmParameter")}: Method
    for generating instances of this class.}
    \item{pData}{\code{signature(object = "HmmParameter")}: See \code{\link[Biobase]{phenoData}}}
    \item{pi}{\code{signature(object = "HmmParameter")}: Accessor fo the
    log initial state probabilities.}
    \item{sampleNames}{\code{signature(object = "HmmParameter")}:
      Accessor for the sample names.}
    \item{show}{\code{signature(object = "HmmParameter")}: Method for
      summarizing the data in instance of this class}
    \item{states}{\code{signature(object = "HmmParameter")}: The hidden
      states.  With the exception of 'N' for normal, there is no
      controlled vocabulary.}
    \item{tau}{\code{signature(object = "HmmParameter")}: Accessor for
      the transition probabilities}
    \item{tau.scale}{\code{signature(object = "HmmParameter")}: See also \code{\link{tau.scale}}}    
	 }
}
\references{

  RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

}
\author{R. Scharpf}
\seealso{
  \code{\link{HmmOptions-class},
    \code{\link[oligoClasses]{oligoSnpSet-class}},
    \code{\link[oligoClasses]{SnpCopyNumberSet-class}},
    \code{\link[oligoClasses]{SnpCallSet-class}}
  }
  }

\examples{
showClass("HmmParameter")
}  
\keyword{classes}
