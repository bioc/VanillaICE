\name{hmm-methods}
\docType{methods}
\alias{hmm}
\alias{hmm-methods}
\alias{hmm,BeadStudioSet-method}
\alias{hmm,CNSet-method}
\alias{hmm,oligoSnpSet-method}
\alias{hmm,SnpSet-method}
\title{Hidden Markov Model methods}
\description{
 Hidden Markov Model methods in package \pkg{VanillaICE}
}
\section{Methods}{


  \describe{

    The \code{hmm} method is defined for several classes of containers
    of preprocessed and normalized SNP array data.  The most common
    containers for use with genotyping platforms are the
    \code{BeadStudioSet} and \code{oligoSnpSet} classes.  The primary
    difference between these two containers are the requirements for the
    assay data elements.  A \code{BeadStudioSet} object must have assay
    data elements "lrr" (log R ratios) and "baf" (B allele frequencies).
    Genotype calls are optional for the \code{BeadStudioSet} object.  As
    the name implies, the \code{BeadStudioSet} container would typically
    be generated as part of a pipeline to process data from Illumina
    array platforms.  By contrast, the \code{oligoSnpSet} object has
    required assay data elements "call" (genotype calls),
    "callProbability" (genotype confidence scores), "copyNumber", and
    "cnConfidence".  As B allele frequencies are perhaps more
    informative than the genotype calls for distinguishing copy number
    states (particularly amplifications), an assay data element named
    "baf" can be included in the assay data for an \code{oligoSnpSet}
    object.  The presence of a "baf" element in the assay data of an
    \code{oligoSnpSet} has implications on the particular HMM fit to
    identify the CNV boundaries (as discussed below).

    \item{\code{signature(object = "BeadStudioSet", ...)}}{

      A hidden Markov model for the \code{BeadStudioSet} class.  The
  assay data are log R ratios and B allele frequencies. See
  \code{hmmBeadStudioSet} for additional arguments that can be passed
  through the \code{...} operator.

}

\item{\code{signature(object = "SnpSet", ...)}}{ A hidden Markov model
  for the \code{SnpSet} class.  The assay data are diallelic genotype
  calls represented as integers (1=AA, 2=AB, 3=BB).  See
  \code{hmmSnpSet} for additional arguments that can be passed
  through the \code{...} operator.
}

\item{\code{signature(object = "CNSet", ...)}}{

  A hidden Markov model for the \code{CNSet} class.  The \code{CNSet}
  instance is first coerced to an object of class \code{oligoSnpSet}
  containing estimates of total copy number and B allele frequencies.
  See \code{hmmBeadStudioSet} for additional arguments that can be
  passed through the \code{...} operator.  For large data sets, the
  initial coercion to the \code{oligoSnpSet} class can be very expensive
  in terms of I/O and require a large amount of RAM.  Users with large
  data sets may prefer to coerce selected samples (e.g., the set of
  samples belonging to a given batch) to an \code{oligoSnpSet} object,
  and then fit the hmm on the \code{oligoSnpSet} object directly.  This
  approach is illustrated in the \code{crlmmDownstream} vignette.

}

\item{\code{signature(object = "CopyNumberSet", ...)}}{

  A hidden Markov model for the \code{CopyNumberSet} class.  The assay
  data are estimates of total copy number.  This method should not be
  used for arrays with genotype information as the genotypes / B allele
  frequencies are informative for copy number inference.

}

\item{\code{signature(object = "oligoSnpSet", ...)}}{

A hidden Markov model for the \code{oligoSnpSet} class.  If "baf" is
included among the assay data elements, the \code{hmmBeadStudioSet} HMM
is implemented.  Otherwise, the \code{hmmOligoSnpSet} is implemented.

}

}
}

\seealso{\code{\link{hmmBeadStudioSet}}, \code{\link{hmmOligoSnpSet}},
  \code{\link{hmmSnpSet}}}
\keyword{methods}

