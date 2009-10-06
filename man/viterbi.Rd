\name{viterbi}
\alias{viterbi}
\title{viterbi algorithm}

\description{

  The Viterbi algorithm for computing the most likely state sequence
  given a model

}
\usage{
viterbi(emission, tau, arm, initialStateProbs, verbose =
FALSE, chromosome, position, sampleNames, locusNames, normalIndex,
returnLikelihood = FALSE, normal2altered=1, altered2normal=1,
altered2altered=1, verbose=TRUE)
}

\arguments{
    \item{emission}{matrix of log emission probabilities (one sample is
      a matrix)}
    \item{tau}{transition probabilities (original scale)}
    \item{arm}{numeric or character string indicating chromosomal arm}
    \item{initialStateProbs}{initial state probabilities (log scale)}    
    \item{verbose}{Logical}
    \item{chromosome}{chromosome}
    \item{position}{physical position}
    \item{sampleNames}{sample labels}
    \item{locusNames}{labels for loci}
    \item{normalIndex}{index corresponding to the normal state. See details}
    \item{returnLikelihood}{whether to return the 'loglikelihood'}
    \item{normal2altered}{factor for scaling the probability of transitioning
      from the normal state to an altered state}
    \item{altered2normal}{factor for scaling the probability of transitioning
      from an altered state to a normal state}
    \item{altered2altered}{factor for scaling the probability of transitioning
      from an altered state to a different altered state}
    \item{verbose}{Logical.  Whether to display all messages and warnings.}
}
\details{

  The Viterbi algorithm is fit independently to each chromosomal arm if
  arm is specified.

}
\value{
  \item{matrix}{predicted states}
}
\author{R. Scharpf}
\keyword{models}

