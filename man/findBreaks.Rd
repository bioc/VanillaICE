\name{findBreaks}
\alias{findBreaks}
\title{Identify breakpoints from the hidden Markov model predictions}
\description{
Identify breakpoints: physical position of breaks, number of SNPs in
region, and the called hidden state.
}
\usage{
findBreaks(x, states, position, chromosome, sample, lik1, lik2)
}
\arguments{
  \item{x}{matrix of predicted hidden states (integer)}
  \item{states}{character string giving the hidden states}
  \item{position}{physical position of the SNPs}
  \item{chromosome}{chromosome number}
  \item{sample}{sample name}
  \item{lik1}{likelihood under alternative model}
  \item{lik2}{likelihood under NULL model}
}
\details{
}
\value{
  \item{data.frame}{data.frame describing regions of normal and altered states}
}
\examples{
}
\keyword{manip}

