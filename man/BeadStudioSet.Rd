\name{BeadStudioSet}
\alias{BeadStudioSet}

\title{
Constructor for BeadStudioSet class
}
\description{
  Constructs an instance of BeadStudioSet from a list of files
  containing log R ratios and B allele frequencies.
}
\usage{
BeadStudioSet(filenames, lrr.colname = "Log.R.Ratio", baf.colname = "B.Allele", sep = "\t", header = TRUE, colClasses, genome = c("hg19", "hg18"), annotationPkg, chromosome = 1:22, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{filenames}{ character string providing the names of the
    BeadStudio files, including the complete path if not in the working
    directory.  }

  \item{lrr.colname}{
    character string providing the column header for log R ratios
  }

  \item{baf.colname}{
    character string providing the column header for log R ratios
  }

  \item{sep}{
    field delimiter in the BeadStudio files. See \code{read.table}
  }

  \item{header}{
    logical: whether the files contain a header.
  }

  \item{colClasses}{
    See \code{read.table}.
  }

  \item{genome}{
    Character string indicating which genome build to use.  Supported entries are "hg19" and "hg18".
}

  \item{annotationPkg}{
    character string providing the name of the annotation package.
  }

  \item{chromosome}{ integer vector indicating which chromosomes to
    include in the BeadStudioSet.  E.g., 1:23 for autosomes and
    chromosome X }

  \item{\dots}{
    Additional arguments to \code{read.bsfiles}.
}

}

\value{
  An object of class \code{BeadStudioSet}
}

\author{
R. Scharpf
}

\seealso{
  \code{\link{read.bsfiles}}, \code{\linkS4class{BeadStudioSet}}
}
\examples{
	path <- system.file("extdata", package="VanillaICE")
	fname <- file.path(path, "LRRandBAF.txt")
	bsSet <- BeadStudioSet(fname, annotationPkg="genomewidesnp6Crlmm", genome="hg19")
}
\keyword{IO}
\keyword{classes}
