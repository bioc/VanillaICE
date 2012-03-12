\name{hmmSnpSet}
\alias{hmmSnpSet}

\title{
  Function for fitting a HMM to SnpSet containers
}

\description{

  Function for fitting a HMM to SnpSet containers. This HMM uses only
  the genotypes to find regions of homozygosity.  For copy number
  inference, see \code{hmmBeadStudioSet} and \code{hmmOligoSnpSet}.

}

\usage{
hmmSnpSet(object, ICE = FALSE, chromosome = 1:22, normalIndex = 1L, rohIndex = normalIndex + 1L, S = 2L, prGtHom = c(0.95, 0.05), prGtMis = rep(1/S, S), prHetCalledHom = 0.001, prHetCalledHet = 0.995, prHomInNormal = 0.8, prHomInRoh = 0.999, TAUP=1e8,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{SnpSet}.
}
  \item{ICE}{
    Whether to use the genotype confidence scores when estimating the
    emission probabilities.
  }

  \item{chromosome}{ Numeric vector indicating which chromosomes to fit
    for the HMM.  See \code{unique(chromosome(object))} for valid
    chromosomes.}

  \item{normalIndex}{
    Index for state with typical rate of heterozygosity.
  }

  \item{rohIndex}{
    Index for state with homozygous genotypes.
  }

  \item{S}{
    Integer indicating number of states (typically 2).
  }

  \item{prGtHom}{
    Numeric vector indicating the probability of a homozygous genotype
    for each of the hidden states.  E.g., c(0.70, 0.99) for states
    corresponding to typical heterozygosity and homozygosity.
  }

  \item{prGtMis}{
    Numeric vector indicating the probability of a missing genotype for
    each hidden state.  The default assumes that missing genotypes are
    equally probable in any of the hidden states.
  }

  \item{prHetCalledHom}{ Numeric vector indicating the probability that
    a true heteroygous genotype is incorrectly called homozygous -- one
    value for each hidden state.}

  \item{prHetCalledHet}{
    Numeric vector indicating the probability that a truly heterozygous
    genotype is correctly called heterozygous -- one value for each
    hidden state.
}

  \item{prHomInNormal}{ The probability of a homozygous genotype in a
    region with typical heterozygosity.  }

  \item{prHomInRoh}{ The probability of a homozygous genotype in a
    region of homozygosity.  }

  \item{TAUP}{ scalar for defining transition probabilities. Larger
  values of TAUP discourage jumps between states.  }

  \item{\dots}{
    Presently ignored
  }
}

\value{
A \code{RangedData}-derived class.
}

\author{
R. Scharpf
}

\seealso{
  \code{\link{hmm}}, \code{\link{hmmBeadStudioSet}}, \code{\link{hmmOligoSnpSet}}
}
\keyword{smooth}

