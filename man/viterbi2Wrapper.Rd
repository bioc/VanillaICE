\name{viterbi2Wrapper}
\alias{viterbi2Wrapper}
\title{
  Wrapper function for fitting the viterbi algorithm
}
\description{
  The viterbi algorithm, implemented in C,  estimates the optimal state
  path as well as the forward and backward variables that are used for
  updating the mean and variances in a copy number HMM.
}
\usage{
viterbi2Wrapper(r, b, pos, is.snp, cnStates, chrom, prOutlierBAF =
list(initial=1e-5, max=1e-3, maxROH=1e-5), p.hom = 0.05, TAUP = 1e+08,
is.log, center = TRUE, limits, initialProb = rep(1/length(cnStates),
length(cnStates)), normalIndex = 3L, nupdates = 10, tolerance = 5,
computeLLR=TRUE, returnEmission=FALSE, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
    \code{matrix} of copy number estimates.
  }

  \item{b}{
    \code{matrix} of B allele frequencies
  }


  \item{pos}{
    \code{integer} vector of genomic position along a chromosome.
  }

  \item{is.snp}{
    \code{indicator} for whether the marker is polymorphic. Must be the
    same length as the number of rows in \code{r} and \code{b}, and the
    same length as the vector \code{pos}.
  }

  \item{cnStates}{
    \code{numeric} vector for the initial copy number state means.
  }

  \item{chrom}{
    \code{integer}: the chromosome.
  }

  \item{prOutlierBAF}{
    A list with elements 'initial', 'max', and 'maxROH' corresponding to
    the initial estimate of the probability that a B allele frequency
    (BAF) is an outlier, the maximum value for this parameter over
    states that do not involve homozygous genotypes, and the maximum
    value over states that assume homozygous genotypes.  This parameter
    is experimental and could be used to fine tune the HMM for different
    platforms. For example, the BAFs for the Affy platform are typically
    more noisey than the BAFs for Illumina.  One may want to set small
    values of these parameters for Illumina (e.g, 1e-5, 1e-3, and 1e-5)
    and larger values for Affy (e.g., 1e-3, 0.01, 1e-3).

  }

\item{p.hom}{

    \code{numeric}: weight for observing homozygous genotypes.  For
    value \code{0}, homozygous genotypes / B allele frequencies have the
    same emission probability in the 'normal' state as in the states
    hemizygous deletion and in copy-neutral region of homozygosity.
    Regions of homozygosity are common in normal genomes. For small
    values of \code{p.hom}, hemizygous deletions will only be called if
    the copy number estimates show evidence of a decrease from normal.

  }

  \item{TAUP}{ \code{numeric}: scalar for the transition probability
    matrix.  Larger values discourage transitions from the normal state.
  }

  \item{is.log}{
    \code{logical}: Whether the copy number estimates in the \code{r}
    matrix are on the log-scale.
}


\item{center}{

    \code{logical}: If TRUE, the copy number estimates for a chromosomal
    arm are recentered such that the median value is the value specified
    for the mean of the normal copy number state.

}


  \item{limits}{

    \code{numeric} vector of length two specifying the range of the copy
    number estimates in \code{r}. Values of \code{r} outside of this
    range are truncated. See \code{copyNumberLimits}.

  }


  \item{initialProb}{

    \code{numeric} vector indicating the initial state probabilities for
    the hidden Markov model.  The length of \code{initialProb} must be
    the same as the length of \code{cnStates}.

}

  \item{normalIndex}{
    \code{integer} specifying the index for the normal state.  Note
    that states must be ordered by the mean of the copy number
    state. E.g., state 1 is homozygous deletion (0 copies), state 2 is hemizygous
    deletion (1 copy), normal (2 copies), ...  In a 6-state HMM,
    normalIndex should be 3.
  }


  \item{nupdates}{

    \code{integer} specifying the maximum number of iterations for
    reestimating the mean and variance for each of the copy number
    states.  The number of iterations may be fewer than \code{nupdates}
    if the difference in the log-likelihood between successive
    iterations is less than \code{tolerance}.
  }


  \item{tolerance}{ \code{numeric} value for indicating convergence of
    the log-likelihood. If the difference in the log-likelihood of the
    observed data given the HMM model at iteration i and i-1 is less
    than tolerance, no additional updates of model parameters using the
    EM algorithm is needed.
  }


%  \item{constraint}{ Named list with 3 elements: 'mus.cn', 'mus.baf',
%    'sds.baf', corresponding to constraints of the copy number means,
%    mean B allele frequency, and B allele frequency standard deviations.
%    See VanillaICE:::constrainMu, VanillaICE:::constrainMuBaf, and
%    VanillaICE:::constrainSdBaf for the default constraints.
%    Experimental.
%
%  }

  \item{computeLLR}{

    Logical.  Whether to compute a log likelihood ratio (LLR) comparing
    the predicted state to the normal state. This is calculated post-hoc
    and is not precisely the likelihood estimated from the Viterbi
    algorithm. When \code{FALSE}, the LLR is not calculated and the
    algorithm is slightly faster.

  }

  \item{returnEmission}{

    Logical. If TRUE, an array of emission probabilities are returned.
    The dimensions of the array are SNPs, samples, and copy number
    states.

  }

  \item{verbose}{

    Logical. Whether to print some of the details of the processing.

    }

  \item{\dots}{
    Additional arguments can be passed to the function
    \code{cnEmissionFromMatrix} and is currently for internal use only.
  }
}

\details{
  This function is used by related packages extending \pkg{VanillaICE}
  and is not intended to be called directly by the user.
}

\value{

  A \code{GRanges} object if \code{returnEmission} is FALSE.  Otherwise,
  an array is returned.

}

\author{
R. Scharpf
}
\keyword{manip}

