\name{HmmOptions-class}
\docType{class}
\alias{HmmOptions-class}
%\alias{.gtEmission.ICE,integer,HmmOptions-method}
\alias{SnpClass,HmmOptions-method}
\alias{annotation,HmmOptions-method}
\alias{annotation<-,HmmOptions,ANY-method}
\alias{calculateTransitionProbability,SnpLevelSet,HmmOptions-method}
\alias{cn.ICE,HmmOptions-method}
\alias{cn.ICE<-,HmmOptions,logical-method}
\alias{cn.location,HmmOptions-method}
\alias{cn.robustSE}
\alias{cn.robustSE,HmmOptions-method}
\alias{cnEmission,numeric,HmmOptions-method}
\alias{db,HmmOptions-method}
\alias{gt.ICE,HmmOptions-method}
\alias{gt.ICE<-,HmmOptions,logical-method}
\alias{gtEmission,numeric,HmmOptions-method}
\alias{gtEmission,integer,HmmOptions-method}
\alias{initialize,HmmOptions-method}
\alias{notes,HmmOptions-method}
\alias{show,HmmOptions-method}
\alias{states}
\alias{states<-}
\alias{states,HmmOptions-method}
\alias{states<-,HmmOptions,character-method}


\title{Options for the hidden Markov model}
\description{Options for fitting the hidden Markov model are stored in
  objects of this class.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HmmOptions",
  snpset, states, cn.location, cn.SE, cn.robustSE, cn.ICE, gt.gte,
  gte.state, gt.confidence.states, gt.ICE, annotation, beta, notes)}.
}
\section{Slots}{
	 \describe{
    \item{\code{SnpClass}:}{Object of class \code{"character"}
      specifying one of the classes for high throughput SNP data (e.g., \code{\link[oligoClasses]{oligoSnpSet-class}}).}
    \item{\code{states}:}{Character string specifying the hidden states
      of the HMM.  With the exception of 'N' for normal, no controlled
      vocabular is required for the altered states }
    \item{\code{cn.location}:}{Numerical vector for the mean of the
      Gaussian distributions -- specified for each hidden state
      (typically on log scale)}
    \item{\code{cn.robustSE}:}{Numeric:  a robust standard error of the
      copy number distribution -- estimated using quantiles of the
      observed copy number distribution}
    \item{\code{cn.ICE}:}{Logical: whether to use confidence scores for
      the genotype calls in the HMM.}
    \item{\code{gt.ICE}:}{Logical: whether to use confidence scores for
      the copy number estimates (inverse standard errors) in the HMM}
    \item{\code{gt.gte}:}{Numerical vector for the conditional
      probability of the true genotype conditional on the observed
      genotype at a SNP (e.g., P(true genotype = AB | observed genotype
      = AA or BB)).  Note used if gt.ICE is FALSE.}
    \item{\code{gte.state}:}{Character vector for the names of the
      hidden states for the genotype data (e.g., L [loss] and N
      [retention or normal])}
    \item{\code{gt.confidence.states}:}{Character vector for the hidden
      states for genotype.}
    \item{\code{annotation}:}{Chip type (e.g, pd.mapping50k.xba240).
      This is a controlled vocabulary.  }
    \item{\code{notes}:}{Character string}
  }
}
\section{Methods}{
  \describe{
%    \item{.gtEmission.ICE}{\code{signature(x = "integer", options =
%    "HmmOptions")}:... }

    \item{SnpClass}{\code{signature(object = "HmmOptions")}: Accessor
      for the character string specifying the name of the SNP data class
      (e.g., \code{oligoSnpSet})}
    
    \item{annotation}{\code{signature(object = "HmmOptions")}: Accessor
    for the chip type.}
    \item{annotation<-}{\code{signature(object = "HmmOptions", value =
    "ANY")}: Replacement method for specifying chip type.}
    \item{calculateTransitionProbability}{See \code{\link{calculateTransitionProbability}}.}
    \item{cn.ICE}{\code{signature(object = "HmmOptions")}: See also \code{\link{cn.ICE}} }
    \item{cn.ICE<-}{\code{signature(object = "HmmOptions", value =
    "logical")}: See also \code{\link{cn.ICE}}}
    \item{cn.location}{\code{signature(object = "HmmOptions")}: See also
      \code{\link{cn.location}}}
    \item{cn.robustSE}{\code{signature(object="HmmOptions")}: Accessor
    for the robust estimate of the copy number standard deviation.
    Estimated using the autosomal chromosomes.}    
    \item{cnEmission}{\code{signature(x = "numeric", options =
    "HmmOptions")}: See also \code{\link{cnEmission}}}
    \item{db}{\code{signature(object = "HmmOptions")}: Retrieves SNP
    annotation using RSQLite databases.  If varlabels 'chromosome' and
    'position' are present in the featureData, this function is not called.}
    \item{gt.ICE}{\code{signature(object = "HmmOptions")}: See also \code{\link{gt.ICE}}}
    \item{gt.ICE<-}{\code{signature(object = "HmmOptions", value =
    "logical")}: see also \code{\link{gt.ICE}}}
    \item{gtEmission}{\code{signature(x = "numeric", options =
    "HmmOptions")}: See also \code{\link{gtEmission}}}
    \item{gtEmission}{\code{signature(x = "integer", options =
    "HmmOptions")}: Replacement method for genotype emission
    probabilities. See also \code{\link{gtEmission}} }
    \item{initialize}{\code{signature(.Object = "HmmOptions")}: an
    instance of the class can be created using \code{new("HmmOptions", ...)}}
    \item{notes}{\code{signature(object = "HmmOptions")}: character string}
    \item{show}{\code{signature(object = "HmmOptions")}}
    \item{states}{\code{signature(object = "HmmOptions")}: character
    string specifying the hidden states.  'N' is required for normal states.}
	 }
}
\author{R. Scharpf}
\examples{
showClass("HmmOptions")
}  
\keyword{classes}
