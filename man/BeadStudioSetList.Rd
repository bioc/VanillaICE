\name{BeadStudioSetList}
\alias{BeadStudioSetList}

\title{Constructor for BeadStudioSetList class.}
\description{
  Reads processed files containing log R Ratios and B allele frequencies
  and construct a BeadStudioList object.
}
\usage{
BeadStudioSetList(fnames, annotationPkg, genome = c("hg19", "hg18"), outdir = ldPath(), sampleIds, phenoData, byArm=FALSE, centromere, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnames}{
    \code{character} vector containing the complete path to files containing log R ratios and BAFs.
}
  \item{annotationPkg}{
    \code{character} string indicating the name of the annotation package.
}
  \item{genome}{
    \code{character} string indicating which genome build.  Supported
    entries are UCSC builds "hg19" and "hg18".
}
  \item{outdir}{
    \code{character} string indicating where to store ff files for
    storing the log R ratios and B allele frequencies.  Ignored if the
    \code{ff} package is not loaded.
}
  \item{sampleIds}{
    \code{character} vector of sample identifiers.  If missing,
    \code{basename(fnames)} is used.
  }

  \item{phenoData}{
    An \code{AnnotatedDataFrame} containing covariates for the samples.
  }

  \item{byArm}{
       Logical.  Whether each element in the list should be a chromosome
       arm.  If TRUE, centromere must be provided. (experimental)
     }

     \item{centromere}{
       \code{data.frame} containing start and stop coordinates of centromeres.
     }

  \item{\dots}{
    Additional arguments passed to the initialization method for \code{BeadStudioSetList}.
}
}

\value{
  A BeadStudioSetList.
}
\author{
R. Scharpf
}
\seealso{
  \code{\linkS4class{BeadStudioSet}}, \code{\linkS4class{BeadStudioSetList}}
}
\examples{
new("BeadStudioSetList")
}
\keyword{classes}

