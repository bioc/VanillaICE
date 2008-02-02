\name{chromosome1}
\alias{chromosome1}
\docType{data}
\title{Simulated chromosome 1}
\description{
 Simulated genotype calls and copy number estimates for chromosome 1
 (Affymetrix 100k platform).
}
\usage{data(chromosome1)}
\format{
}

\details{

  The simulation comprises one subject's genotype, copy number, and
confidences scores for 9165 SNPs on chromosome 1.  A description of the
5 features simulated in chromosome 1, referred to by regions A-E, and
the underlying hidden states in these regions follows.

Genotype calls: With the exception of three regions in this chromosome,
we simulated 9165 genotypes (the approximate number of SNPs in the two
50k SNP chips) from a Bernoulli distribution with probability 0.7 of a
homozygous genotype.  Confidence scores for genotype estimates were
obtained by random draws of confidence scores in the Hapmap data when
the CRLMM genotype calls agreed with the gold-standard as defined by
consensus of the HapMap genotyping centers. 

Copy number: The Affymetrix CNAT tool (version 3.0) was used to obtain
copy number estimates for the 9165 SNPs from a presumably normal
individual in the HapMap dataset (sample identifier NA06993).  Deletions
and amplificatons were simulated from Gaussian distributions with
location parameters log2(1) and log2(3), respectively.  For the scale
parameter, we used a robust estimate of the log2 transformed copy number
standard deviation.  To illustrate how a confidence score such as a
standard error of the copy number estimate could be useful, we simulated
standard errors from a shifted Gamma: Gamma(1, 2) + 0.3, where 1 is the
shape parameter and 2 is the rate parameter.  To ascertain the effect of
qualitatively high confidence scores on the ICE HMM, we scaled the
robust standard deviation estimate by 1/2. Similarly, to simulate less
precise copy number estimates, we multiplied the robust standard
deviation estimate by a factor of 2.

For additional information on the five abnormalities simulated in this
chromosome, see the manuscript referenced below.
 }

\source{
}

\references{

    RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

}

\examples{
data(chromosome1)
chromosome1
}
\keyword{datasets}
