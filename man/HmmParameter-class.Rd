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
\alias{cn.robustSE,HmmParameter-method}
\alias{fData,HmmParameter-method}
\alias{featureData,HmmParameter-method}
\alias{featureNames,HmmParameter-method}
\alias{gt.ICE,HmmParameter-method}
\alias{hmm,HmmParameter-method}
\alias{hmmOptions}
\alias{hmmOptions<-}
\alias{hmmOptions,HmmParameter-method}
\alias{hmmOptions<-,HmmParameter,HmmOptions-method}
\alias{initialize,HmmParameter-method}
\alias{pData,HmmParameter-method}
\alias{pi}
\alias{pi,HmmParameter-method}
\alias{sampleNames,HmmParameter-method}
\alias{show,HmmParameter-method}
\alias{states,HmmParameter-method}
\alias{states<-,HmmParameter,character-method}
\alias{tau}
\alias{tau,HmmParameter-method}
\alias{tau.scale,HmmParameter-method}

\title{Parameters for fitting a hidden Markov model to SNP chip data}

\description{This class contains all the parameters needed to fit the
  hidden Markov model.}

\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("HmmParameter",
    snpset, states, tau, tau.scale, pi, beta, notes, SCALE, ...)}.
  
}
\section{Slots}{
  \describe{
    \item{\code{hmmOptions}:}{Options for fitting the HMM are stored in an object of class \code{HmmOptions}}
    \item{\code{tau}:}{Transition probabilities. Object of class \code{"matrix"}.  }
    \item{\code{tau.scale}:}{A matrix to scale the transition probabilities.  See also \code{\link{tau.scale}}}
    \item{\code{pi}:}{A numerical vector specifying the log initial state probabilities}
    \item{\code{beta}:}{A matrix of emission probabilities for the observed data.  Each column corresponds to a sample.}
    \item{\code{featureData}:}{Object of class \code{AnnotatedDataFrame}}
    \item{\code{notes}:}{Object of class \code{character}}
  }
}
\section{Methods}{
  \describe{
    \item{Beta}{\code{HmmParameter}: Accessor for  the emission probabilities.  See also \code{\link{Beta}}}
    \item{SnpClass}{\code{HmmParameter}: A character string specifying
      one of the data classes for high throughput SNP data.}
    \item{"["}{\code{HmmParameter}: For subsetting.}
    \item{arm}{\code{HmmParameter}: See also \code{\link{arm}}}
    \item{chromosome}{\code{HmmParameter}: Accessor for chromosome information}
    \item{cn.ICE}{\code{HmmParameter}: See also \code{\link{cn.ICE}}}
    \item{cn.location}{\code{HmmParameter}: See  also \code{\link{cn.location}}}
    \item{cn.robustSE}{\code{HmmParameter}: Accessor for a the robust estimate of the copy number standard deviation   estimated from the autosomal chromosomes}
    \item{fData}{\code{HmmParameter}: Accessor for the data frame of all feature-level (SNP-level) annotation (e.g., physical position)}    
    \item{featureData}{\code{HmmParameter}: Accessor for the \code{\link[Biobase]{featureData}}}
    \item{featureNames}{\code{HmmParameter}:  Accessor for the \code{featureNames}}
    \item{gt.ICE}{\code{HmmParameter}: See also \code{\link{gt.ICE}}}
    \item{hmm}{\code{HmmParameter}: See also \code{\link{hmm}}}
    \item{hmmOptions}{\code{HmmParameter}: Accessor for the options used to fit the HMM.  See also \code{\link{HmmOptions-class}}}
    \item{hmmOptions<-}{\code{HmmParameter, HmmOptions}: Replace the slot hmmOptions in an object of class \code{HmmParameter} with an object of class \code{HmmOptions}}
    \item{initialize}{\code{HmmParameter}: Method  for generating instances of this class.}
    \item{pData}{\code{HmmParameter}: See \code{\link[Biobase]{phenoData}}}
    \item{pi}{\code{HmmParameter}: Accessor fo the log initial state probabilities.}
    \item{sampleNames}{\code{HmmParameter}:  Accessor for the sample names.}
    \item{show}{\code{HmmParameter}: Method for   summarizing the data in instance of this class}
    \item{states}{\code{HmmParameter}: The hidden   states.  With the exception of 'N' for normal, there is no   controlled vocabulary.}
    \item{tau}{\code{HmmParameter}: Accessor for   the transition probabilities}
    \item{tau.scale}{\code{HmmParameter}: See also \code{\link{tau.scale}}}    
  }
}
\references{

  RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

}
\author{R. Scharpf}
\seealso{
  \code{\link{HmmOptions-class}},
    \code{\link[oligoClasses]{oligoSnpSet-class}},
    \code{\link[oligoClasses]{SnpCopyNumberSet-class}},
    \code{\link[oligoClasses]{SnpCallSet-class}}
  }

\examples{
showClass("HmmParameter")
}  
\keyword{classes}
