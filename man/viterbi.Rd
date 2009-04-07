\name{viterbi}
\alias{viterbi}
\title{viterbi algorithm}

\description{

  The Viterbi algorithm for computing the most likely state sequence
  given a model

}
\usage{
viterbi(initialStateProbs, emission, tau, arm, tau.scale, verbose =
FALSE, chromosome, position, sampleNames, locusNames, normalIndex, returnLikelihood = FALSE)
}

\arguments{
    \item{initialStateProbs}{initial state probabilities (log scale)}
    \item{emission}{matrix of log emission probabilities (one sample is a matrix)}
    \item{tau}{transition probabilities (original scale)}
    \item{arm}{numeric or character string indicating chromosomal arm}  
    \item{tau.scale}{matrix to scale the probability of transitioning
      between states.}
    \item{verbose}{Logical}
    \item{chromosome}{chromosome}
    \item{position}{physical position}
    \item{sampleNames}{sample labels}
    \item{locusNames}{labels for loci}
    \item{normalIndex}{index corresponding to the normal state. See details}
    \item{returnLikelihood}{whether to return the 'loglikelihood'}
}
\details{

  The Viterbi algorithm is fit independently to each chromosomal arm if
  arm is specified.

  Argument \code{tau.scale} is a matrix that scales the probability of
  transitioning from an altered state to a normal state to the
  probability of transitioning between two altered states.  If missing,
  tau.scale is 1 (no scaling)

  

}
\value{
  \item{matrix}{predicted states}
}
\author{R. Scharpf}
\keyword{models}

