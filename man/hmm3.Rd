\name{hmm3}
\alias{hmm3}

\title{
  HMM for processed data.
}


\description{

  HMM for preprocessed and normalized data. Expects files containing log
  R ratios and B allele frequencies that have already been
  wave-corrected, etc.

}

\usage{
hmm3(filenames, annotationPkg, universe=c("", "hg18", "hg19"),
lrr.colname = "Log.R.Ratio", baf.colname = "B.Allele", sep="\t",
TAUP=1e8, cnStates = c(-1.5, -0.5, 0, 0, 0.4, 0.8), normalIndex=3L,
rohIndex=normalIndex+1L, prOutlierBAF=1e-3, prOutlierCN=0.01, p.hom = 1,
is.log=TRUE, initialProb=rep(1/length(cnStates), length(cnStates)),
center=TRUE, nupdates=10, tolerance=1, chromosome=1:22, ...)
}

\arguments{
  \item{filenames}{
    character string of filenames:  must include the complete path such
  that \code{file.exists(filenames)} evaluates to \code{TRUE}.
}

\item{annotationPkg}{

  character string identifying annotation package.

}

  \item{universe}{
    Not currently implemented.
}


  \item{lrr.colname}{

    character: used to grep for the log R ratios in the header. E.g.,
    \code{grep(lrr.colname, header) } should return a length 1 vector,
    where header is a vector of the column labels.

  }

  \item{baf.colname}{

    character: used to grep for the B allele frequency in the
    header. E.g., \code{grep(baf.colname, header) } should return a
    length 1 vector, where header is a vector of the column labels.

  }

  \item{sep}{the field delimiter passed to \code{read.table}.}

  \item{TAUP}{Scalar for transition probabilities.  Large values make the probability of transitioning
    between two different states smaller.  The formula for computing the
    transition probability is

    \code{exp (-2*diff(position)/TAUP)},

    where position is the physical position on the genome. The default
    value is 10^8.}

  \item{cnStates}{

    A numeric vector specifying the initial values of the mean copy
    number for each state -- the Gaussian component of the
    Gaussian-uniform mixture. The initial values are the same for each
    sample, but updated so that the estimates are sample- and
    chromosome-specific.  The copy number states are homozygous deletion
    (total CN=0), hemizygous deletion (total CN=1), normal (total CN=2),
    region of homozygosity (total CN=2) single-copy duplication (total
    CN=2), and two-copy duplication (total CN=4). The reestimation
    procedure is described by L. Rabiner, 1989.

  }

  \item{normalIndex}{An integer indicating the state index for the
    normal state. Currently, this is required to be \code{3L}.}

  \item{rohIndex}{An integer indicating the state index for copy-neutral
  homozygosity (copy number 2, no heterozygosity). Currently, this is
  required to be \code{4L}.}

  \item{prOutlierBAF}{

    A numeric value indicating the probability that an observation was
    emitted from the Uniform component of the Gaussian-Uniform mixture
    used to model B allele frequencies. The initial value is assumed to
    be the same for each sample and the same for each copy number state.
    Currently, the mixture probabilities are not updated.

  }

  \item{prOutlierCN}{

    A numeric value indicating the probability that an observation was
    emitted from the Uniform component of the Gaussian-Uniform
    mixture. The initial value is assumed to be the same for each sample
    and the same for each copy number state.  The mixture probabilities
    are updated using the reestimation procedure described previously
    (Rabiner, 1989).  After updating, the mixture probabilities are
    state- and sample-specific.

  }

  \item{p.hom}{

    The weight used to average the probability that the BAF was emitted
    from a region of homozygosity or a region with normal
    heterozygosity.  (Note that the probability of observing a BAF of a
    zero or 1 in a region of homozygosity is slightly higher than the
    probability of observing these values in a normal region. )  Small
    values of \code{p.hom} (e.g., less than 0.5), has the effect of
    making the emission probability for the normal copy number state
    more similar to state '4' for homozygous regions of the genome.  The
    consequence is that longer tracts of homozygosity are required
    before the HMM will transition to state '4'.  As an extreme example,
    setting \code{p.hom = 0} should restrict calls to normal (state 3)
    and altered copy number (states 1-2 and 5-6).

  }

  \item{is.log}{\code{logical} indicating whether the copy number
  estimates in the assay data are log R ratios (normal state has a mean
  near zero) or absolute (normal state has a mean near 2).}

\item{initialProb}{\code{numeric} vector specifying the probability that
  the first marker on a chromosomal arm belongs to each of the possible
  states.  The probabilities must sum to one.}

\item{center}{\code{logical}: whether to center the copy number
  estimates at the mean specified for the normal copy number state
  (e.g., 2 for absolute copy number)}

\item{nupdates}{\code{integer}: maximum number of times the Viterbi
  algorithm is run to update estimates of the mixture probabilities,
  copy number state means, and variances.  If the difference in the log
  likelihood at update i and update i -1 is less than 1 for all samples,
  only i updates is required.}

\item{tolerance}{\code{numeric} The mixture probabilities, copy number
  means, and copy number variances for copy number states 1-6 are
  updated using the reestimation procedure described in
  Rabiner,1989. The maximum number of update is
  \code{nupdates}. However, if the difference in the log likelihood
  comparing the model at update i and the model at update i-1 is less
  than \code{tolerance} for all samples, the number of updates is
  \code{i}. }

\item{chromosome}{\code{integer} vector indicating which chromosomes to process.}
\item{\dots}{

    Currently ignored.

  }
}


\value{

  A object of class \code{RangedDataHMM}.

}

\author{
R. Scharpf
}

\seealso{
\code{\linkS4class{RangedDataHMM}}
}

\examples{
	path <- system.file("extdata", package="VanillaICE")
	filename <- list.files(path, pattern="LRRandBAF", full.names=TRUE)
	rd <- hmm3(filename, annotationPkg="genomewidesnp6Crlmm", universe="")
}
\keyword{smooth}
\keyword{models}
