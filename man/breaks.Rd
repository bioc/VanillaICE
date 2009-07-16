\name{breaks}
\alias{breaks}
\title{Identify breakpoints from the hidden Markov model predictions}
\description{
Identify breakpoints: physical position of breaks, number of SNPs in
region, and the called hidden state.
}
\usage{
breaks(x, states, position, chromosome, chromosomeAnnotation = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Locus X sample matrix of hidden states where the hidden
    states are represented as integers}
  \item{states}{Labels for the hidden states}
  \item{position}{Physical position of loci}
  \item{chromosome}{integer indicating chromosome (23=X)}
  \item{chromosomeAnnotation}{chromosome annotation. see details}
  \item{verbose}{verbose output}
}
\details{
  One may provide their own chromosome annotation with centromere start
  and stop sites.  The format must be the same as the
  chromosomeAnnotation dataset in the R package SNPchip.
}
\value{
  data.frame

  \item{sample}{sample label}
  \item{chr}{chromosome (23 = X)}
  \item{start}{starting physical position of segment} 
  \item{end}{last physical position of segment}
  \item{nbases}{number of bases in segment}
  \item{nprobes}{number of probes in segment}
  \item{state}{label for the state of the segment}
  
}
\author{R. Scharpf}
\examples{
x <- matrix(rep(c(1, 2, 3, 1, 2), each=50), ncol=1)
breaks(x, states=c("0", "1", "2"), position=1:nrow(x), chromosome=1)
}
\keyword{manip}

