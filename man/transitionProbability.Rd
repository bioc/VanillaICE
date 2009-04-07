\name{transitionProbability}
\alias{transitionProbability}
\title{Compute the transition probability}
\description{
Wrapper for computing the locus-specific transition probability
}
\usage{
transitionProbability(chromosome, position, TAUP = 1e+08, chromosomeAnnotation, verbose = FALSE)
}
\arguments{
  \item{chromosome}{chromosome (integer representation)}
  \item{position}{physical position}
  \item{TAUP}{Scalar for computing transition probabilities (see Details).}
  \item{chromosomeAnnotation}{Optional: chromosome annotation}
  \item{verbose}{Logical: verbose output}
}
\details{

The HMM uses locus-specific transition probabilities that are calculated
as a function of the physical distance between loci. Specifically, the
probability that the locus at position \eqn{t-1} is not informative for
the locus at position \eqn{t} is calculated as \eqn{1-exp(-d/TAUP)},
where d is the physical distance between locus t and locus t-1.  The
default for TAUP is \eqn{1 x 10^8} and can be specified to acheive a
desired amount of sensitivity and specificity.  Larger values of TAUP
decreases the probability of transitioning to other states, and
therefore provides a more smooth fit.

  
}
\value{

  The transitionProbability function (i) transforms the physical
distance between adjacent loci to an estimate of the genomic distance
and (ii) adds an 'arm' variable to the annotation matrix.

\item{chromosome}{chromosome }
\item{position}{physical position}
\item{arm}{an integer.  The HMM uses the arm variable as a factor and is
fit independently to each 'arm'.}
\item{transitionPr}{ locus-specific transition probabilities}
}
\author{R. Scharpf}
\seealso{\code{\link[SNPchip]{chromosomeAnnotation}}}
\keyword{arith}
