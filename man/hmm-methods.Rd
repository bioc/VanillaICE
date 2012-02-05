\name{hmm-methods}
\docType{methods}
\alias{hmm}
\alias{hmm-methods}
%\alias{hmm,BeadStudioSet,ANY-method}
%\alias{hmm,CNSet,HmmOptionList-method}
\alias{hmm,BeadStudioSet-method}
\alias{hmm,CNSet-method}
%\alias{hmm,CopyNumberSet,HmmOptionList-method}
%\alias{hmm,oligoSnpSet,ANY-method}
\alias{hmm,oligoSnpSet-method}
%\alias{hmm,SnpSet,HmmOptionList-method}
\alias{hmm,SnpSet-method}
%\alias{hmm,SnpSet,ANY-method}
\title{Hidden Markov Model methods}
\description{
 Hidden Markov Model methods in package \pkg{VanillaICE}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "BeadStudioSet", hmm.params = "ANY")}}{
  A hidden Markov model for the \code{BeadStudioSet} class.  The assay
  data are log R ratios and B allele frequencies.
}

\item{\code{signature(object = "SnpSet", hmm.params = "ANY")}}{ A hidden
  Markov model for the \code{SnpSet} class.  The assay data are
  diallelic genotype calls represented as integers (1=AA, 2=AB, 3=BB).
  }

\item{\code{signature(object = "CNSet", hmm.params = "HmmOptionList")}}{
  A hidden Markov model for the \code{CNSet} class.  The \code{CNSet}
  instance is first coerced to an object of class \code{oligoSnpSet}
  containing estimates of total copy number and B allele frequencies.
  See the \code{hmm} method for \code{oligoSnpSet} objects.
}

\item{\code{signature(object = "CopyNumberSet", hmm.params = "HmmOptionList")}}{
  A hidden Markov model for the \code{CopyNumberSet} class.  The assay
  data are estimates of total copy number.  This method should not be
  used for arrays with genotype information as the genotypes / B allele
  frequencies are informative for copy number inference.
}

\item{\code{signature(object = "oligoSnpSet", hmm.params = "ANY")}}{

A hidden Markov model for the \code{CopyNumberSet} class.  The assay
  data are estimates of total copy number and genotypes.  If B allele
  frequencies are available (preferred), see the \code{hmm} method for
  code{BeadStudioSet} objects.

}

\item{\code{signature(object = "SnpSet", hmm.params = "HmmOptionList")}}{

  A hidden Markov model for the \code{SnpSet} class.  The assay data are
  estimates of the genotypes.  The hidden states are regions of
  homozygosity and normal heterozgosity.

}

}}
\keyword{methods}

