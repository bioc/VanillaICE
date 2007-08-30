\name{distance}
\alias{distance}
\alias{distance<-}

\title{Accessor for distance}

\description{ Indicates whether distance between SNPs (in base pairs)
was used to estimate transition probabilities for the HMM}

\usage{
distance(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ one of the HMM classes defined in \code{ICE}}
}
\details{
}
\value{
  Logical: if TRUE, the distance between SNPs was used to derive a
  SNP-specific transition probability matrix.

}

\references{

  \cite{RB Scharpf et al. (2007), Department of Biostatistics
    Working Papers, Working Paper 136}

  \cite{MA Newton et al. (1998) On the statistical analysis of
  allelic-loss data. Stat Med 17(13):1425-45}

  }

\author{}
\note{ 


}
\seealso{}
\examples{
}
\keyword{methods}

