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
viterbi2Wrapper(r, b, gt, pos, is.snp, cnStates, chrom, prOutlierBAF = 0.001, p.hom = 0.05, TAUP = 1e+08, is.log, center = TRUE, reestimation = TRUE, limits, initialProb = rep(1/length(cnStates), length(cnStates)), normalIndex = 3L, rohIndex = normalIndex + 1L, nupdates = 10, tolerance = 1, returnViterbiObject = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
    \code{matrix} of copy number estimates.
}
  \item{b}{
    \code{matrix} of B allele frequencies
}
  \item{gt}{
    \code{matrix} of genotype calls (1=AA, 2=AB, 3=BB).  Ignored unless
    \code{b} is missing.
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
    \code{numeric}: initial probability for observing an outlier in the
    B allele frequencies.
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

  \item{reestimation}{

    \code{logical}: if TRUE, the initial values provided for the mean of
    the copy number states will be reestimated as described previously
    (Rabiner, 1989) as a function of the forward and backward
    probabilities from the Viterbi algorithm and the mixture
    probabilities from the Normal-uniform mixture model for each state.

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
    normalIndex should be 3. }

  \item{rohIndex}{
    \code{integer} specifying the index for copy-neutral region of
    homozygosity.  In a 6-state HMM, the \code{rohIndex} should be 4.
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
    EM algorithm is needed.  }

  \item{returnViterbiObject}{
    \code{logical}: whether to return an object of class
    \code{Viterbi}.  For internal use only.
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
  A \code{RangedDataHMM} object if \code{returnViterbiObject} is FALSE.
}

\author{
R. Scharpf
}
\keyword{manip}

