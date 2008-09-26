\name{genotypeEmission}
\alias{genotypeEmission}
\title{Emission probabilities for di-allelic genotype calls}
\description{
Emission probabilities for di-allelic genotype calls
}
\usage{
genotypeEmission(genotypes, states, probHomCall, probMissing, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{matrix of integers (1=AA, 2=AB, 3=BB, 4=other)}
  \item{states}{character string of hidden states}
  \item{probHomCall}{numeric: probability of a homozygous genotype call
    specified in the same order as the hidden states}
  \item{probMissing}{numeric: probability of a missing genotype call
    specified in the same order as the hidden states}
  \item{verbose}{logical}
}
\details{
  By default, missing genotype calls will be assumed to be independent
  of the hidden state.  This is not always appropriate, but is easiest
  to implement in terms of a default. 
}
\value{
  \item{array}{Array of emission probabilities.  Dimension 1: SNPs, Dimension 2: samples, Dimension3: states}
}
\keyword{methods}

