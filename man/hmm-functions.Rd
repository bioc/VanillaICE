\name{hmm-functions}
\alias{hmmBafLrrSet2}
\alias{hmmOligoSnpSet2}
\alias{hmmBafLrrSetList2}

%\alias{hmmOligoSetList2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   HMM functions for oligoSnpSet and BeadStudioSet containers
}
\description{

   HMM functions for \code{oligoSnpSet}, \code{BafLrrSet},
   \code{SnpSet2} and \code{BafLrrSetList} containers. These functions
   are exported in the package's namespace to provide documentation of
   arguments that can be passed from the \code{hmm} method for these
   containers.  The \code{hmmBeadStudioSet2} function is always called
   when the object passed to the \code{hmm} method is a
   \code{BeadStudioSet} or \code{BafLrrSet}.  The \code{hmm} method for
   \code{oligoSnpSet} objects will also call the \code{hmmBafLrrSet2}
   function if if B allele frequencies (assay data element "baf") is
   included in the list of assay data elements. Otherwise, the
   \code{hmm} method for \code{oligoSnpSet} will call the
   \code{hmmOligoSnpSet2} function.

}
\usage{

hmmBafLrrSet2(object, sampleIds, TAUP=1e10, tauMAX,cnStates=c(-2, -0.4, 0, 0, 0.4, 1), is.log=TRUE,...)
hmmOligoSnpSet2(object, sampleIds, TAUP=1e10, tauMAX, cnStates = c(0,1,2,2,3,4), is.log=FALSE, ...)
hmmBafLrrSetList2(object, sampleIds, TAUP=1e10, tauMAX, cnStates=c(-2, -0.4, 0, 0, 0.4, 1), is.log=TRUE,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    Any of the containers listed above.
  }

  \item{sampleIds}{A character vector indicating which samples in
  \code{object} to process. For large data, specifying \code{sampleIds} is
  much more efficient that subsetting the \code{object} prior to fitting
  the HMM.}

\item{TAUP}{ Scalar for the probability that the state of segment t is
  the same as the state at segment t-1.  Larger values decrease the
  probability of transitioning to an different state.}

\item{tauMAX}{Upper bound for the probability that the state at marker t
  is the same as the state at marker t-1.}

\item{cnStates}{

  A vector of starting values (numeric) specifying the means of the
  Normal distribution assumed for latent copy numbers. The means must be
  specified for states homozygous deletion (zero copies), hemizygous
  deletion (1 copy), normal (2 copies), normal and no heterozygotes (2
  copies), single copy duplication (3 copies), and two+ copy duplication
  (4+ copies). The starting values are updated via EM.

}

%\item{normalIndex}{
%
%  Integer indicating the state index for diploid copy number. This
%    should nearly always be '3' if the 6-state HMM (see \code{cnStates})
%    is fit as recommended.
%
%  }


%\item{prOutlierCN}{
%    The probability that a copy number estimate is an outlier.  This is
%    an initial estimate that is updated for each copy number state via
%    EM.
%}

%  \item{p.hom}{
%    \code{numeric}: weight for observing homozygous genotypes.  For
%    value \code{0}, homozygous genotypes / B allele frequencies have the
%    same emission probability in the 'normal' state as in the states
%    hemizygous deletion and in copy-neutral region of homozygosity.
%    Regions of homozygosity are common in normal genomes. For small
%    values of \code{p.hom}, hemizygous deletions will only be called if
%    the copy number estimates show evidence of a decrease from normal.
%
%  }
%  \item{TAUP}{
%
%    Scalar for the transition probability matrix.  Larger values
%    discourage transitions from the normal state. (The transition
%    probabilities are a function of the distance between adjacent
%    markers.  These probabilities are not updated as part of the EM
%    step.)
%
%  }

  \item{is.log}{

    A \code{logical} indicating whether the copy number estimates are on
    the log scale.  Note that the assay data elements in
    \code{oligoSnpSet} and \code{BeadStudioSet} should be represented as
    integers (copy number or relative copy number * 100). If
    \code{is.log} is TRUE, we assume that after division by 100 the
    assay data element containing the copy numbers (or relative copy
    numbers) is on the log-scale.  The scale has implications on what is
    considered to be extreme.

  }

%  \item{initialProb}{
%
%    Vector of initial state probabilities.  This is required to be the
%    same length as \code{cnStates}.
%
%  }

%  \item{center}{
%
%    Whether to center the copy number for each chromosomal arm at the
%    theoretical mean for the diploid copy number state.  This may not be
%    appropriate for some datasets (e.g., trisomy 21, cancer
%    applications).  A safer approach is to set this argument to
%    \code{FALSE} and center all autosomes at the theoretical mean for
%    two copies prior to fitting the HMM.
%
%  }

%  \item{nupdates}{
%
%    The maximum number of reestimation steps for updating the mean,
%    variance, and outlier probabilities of the Gaussian-Uniform mixture
%    for each copy number state.
%
%  }

%  \item{tolerance}{
%
%    If the difference in the log likelihood between successive EM
%    updates is less than tolerance, the number of updates can be less
%    than \code{nupdates}.
%
%  }

%  \item{sampleIds}{ Character vector indicating which samples to
%    process.  If missing, all samples in \code{object} are processed.}
%
%  \item{computeLLR}{
%
%    Logical.  Whether to compute a log likelihood ratio (LLR) comparing
%    the predicted state to the normal state. This is calculated post-hoc
%    and is not precisely the likelihood estimated from the Viterbi
%    algorithm. When \code{FALSE}, the LLR is not calculated and the
%    algorithm is slightly faster.
%
%  }
%

  \item{\dots}{
    Additional arguments can be passed to \code{viterbi2Wrapper}.
  }
}
\value{
  A \code{GRanges} or \code{GRangesList} object of the segmented data.
}

\author{
R. Scharpf
}

%\note{ \code{prOutlierBAF} can be passed to the \code{viterbi2Wrapper}
%  function through the \code{...} operator. It may be desirable to
%  specify different values for this parameter depending on whether the
%  platform is Affymetrix or Illumina.  See \code{viterbi2Wrapper} for
%  additional details.  }

%\seealso{
%  \code{\link{hmmSnpSet}}, \code{\link{viterbi2Wrapper}},
%  \code{\link{hmm}}, \code{\link{hmm-methods}}
%}
\seealso{hmm-methods}

\keyword{smooth}

