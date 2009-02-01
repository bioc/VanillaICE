\name{genotypeEmission}
\alias{genotypeEmission}
\title{Emission probabilities for di-allelic genotype calls}
\description{
Emission probabilities for di-allelic genotype calls
}
\usage{
genotypeEmission(genotypes, conf, states, probHomCall, probMissing, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{matrix of integers (1=AA, 2=AB, 3=BB, 4=other)}
  \item{conf}{Confidence estimates of the genotype calls obtained from crlmm (optional).}
  \item{states}{character string of hidden states}
  \item{probHomCall}{numeric: probability of a homozygous genotype call
    specified in the same order as the hidden states}
  \item{probMissing}{numeric: probability of a missing genotype call
    specified in the same order as the hidden states}
  \item{verbose}{logical}
}

\details{
  
  CRLMM provides confidences estimates of the genotype calls that can be
  integrated to improve the HMM.  Because CRLMM will genotype all SNPs,
  the probMissing argument is unecessary.
  
}


\value{
  \item{array}{Array of emission probabilities.  Dimension 1: SNPs, Dimension 2: samples, Dimension3: states}
}

\keyword{methods}

