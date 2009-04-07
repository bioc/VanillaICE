\name{hmm}
\alias{hmm}
\title{Wrapper for fitting the HMM}
\description{
  A wrapper for fitting the HMM.
}
\usage{
hmm(object, states, mu = NULL, probs = NULL, takeLog = FALSE, initialP, returnSegments = TRUE, TAUP = 1e+08, verbose = FALSE, ice = FALSE, envir)
}

\arguments{
    \item{object}{SnpCallSet, SnpCopyNumberSet, or oligoSnpSet object}
  \item{states}{Labels for the hidden states.  See details for order.}
  \item{mu}{The latent copy number. See details for order. }
  \item{probs}{See details.}
  \item{takeLog}{Whether to take the log of the copy number before
    computing emission probabilities and standard deviations}
  \item{initialP}{Initial state probabilities}
  \item{returnSegments}{Logical: whether to return the segments or the
    loci x sample matrix of predicted states}
  \item{TAUP}{Scaling parameter for transition probabilities.}
  \item{verbose}{Logical: Verbose output?}
  \item{ice}{Whether to use CRLMM confidence scores of the genotype calls.}
  \item{envir}{Optional.  An environment for storing intermediate files
    created for fitting the HMM.}
}
\details{
  
  For oligoSnpSet objects, the hidden state labels are assumed to be
  1: hemizygous deletion
  2: normal
  3: region of homozygosity (ROH)
  4: amplification

  The argument mu should have copy number values corresponding to the
  above states.  For instance on the absolute scale, the copy number
  states should be 1, 2, 2, and 4.

  \code{probs}:  If \code{ice} is FALSE, the elements in probs should
  correspond to the probability of a homozygous genotype in each of the
  above states.  If \code{ice} is TRUE, the elements in probs should
  correspond to
  1.  Pr(homozygous call | truth is heterozyous)  
  2.  Pr(heterozygous call | truth is heterozygous)
  3.  Pr(homozygous call | truth is ROH)
  4.  Pr(homozygous call | truth is normal) .  'Normal' meaning copy
  number 2 and a typical frequency of heterozygosity for autosomes.
}
\value{

  If \code{returnSegments} is TRUE, a data.frame containing the
  coordinates of the predicted segments is returned.  Otherwise, a loci
  X sample matrix is returned. The elements of the matrix correspond to
  the predict hidden state for a specific locus and sample.
}

\references{

    RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

  }
\author{R. Scharpf}
\keyword{models}
