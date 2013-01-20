\name{hmmSnpSet}
\alias{hmmSnpSet2}
\alias{hmmSnpSetIce}
\title{
  Function for fitting a HMM to SnpSet containers
}

\description{

  Function for fitting a HMM to SnpSet containers. This HMM uses only
  the genotypes to find regions of homozygosity.  For copy number
  inference, see \code{hmmBeadStudioSet} and \code{hmmOligoSnpSet}.

}

\usage{
hmmSnpSet2(object, sampleIds, TAUP=1e10, tauMAX,normalIndex=1L,rohIndex=normalIndex+1L,S=2L,...)
hmmSnpSetIce(object, sampleIds, TAUP=1e10, tauMAX, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{SnpSet} or \code{SnpSet2} object.
  }

  \item{sampleIds}{A character vector indicating which samples in
  \code{object} to process. For large data, specifying \code{sampleIds} is
  much more efficient that subsetting the \code{object} prior to fitting
  the HMM.}

  \item{TAUP}{ scalar for defining transition probabilities. Larger
    values of TAUP discourage jumps between states.  }

  \item{tauMAX}{Upper bound for the probability that the state at marker t
  is the same as the state at marker t-1.}

  \item{normalIndex}{
    Index for state with typical rate of heterozygosity.
  }

  \item{rohIndex}{
    Index for state with homozygous genotypes.
  }

  \item{S}{
    Integer indicating number of states (typically 2).
  }

%  \item{prGtHom}{
%    Numeric vector indicating the probability of a homozygous genotype
%    for each of the hidden states.  E.g., c(0.70, 0.99) for states
%    corresponding to typical heterozygosity and homozygosity.
%  }
%
%  \item{prGtMis}{
%    Numeric vector indicating the probability of a missing genotype for
%    each hidden state.  The default assumes that missing genotypes are
%    equally probable in any of the hidden states.
%  }
%
%  \item{prHetCalledHom}{ Numeric vector indicating the probability that
%    a true heteroygous genotype is incorrectly called homozygous -- one
%    value for each hidden state.}
%
%  \item{prHetCalledHet}{
%    Numeric vector indicating the probability that a truly heterozygous
%    genotype is correctly called heterozygous -- one value for each
%    hidden state.
%}
%
%  \item{prHomInNormal}{ The probability of a homozygous genotype in a
%    region with typical heterozygosity.  }
%
%  \item{prHomInRoh}{ The probability of a homozygous genotype in a
%    region of homozygosity.  }



  \item{\dots}{
    Additional arguments passed to \code{viterbiForSnpSetIce}
  }
}

\value{
A \code{GRanges} or \code{GRangesList} object containing the genomic
intervals from the HMM annotated with the copy number states and number
of probes per interval.
}

\author{
R. Scharpf
}

\seealso{
  \code{\link{hmm}}, \code{\link{hmmBafLrrSet2}}
}
\keyword{smooth}

