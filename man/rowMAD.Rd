\name{rowMAD}
\alias{rowMAD}

\title{
  Calculate the median absolute deviation for each row in a matrix.
}
\description{
  Calculate the median absolute deviation for each row in a matrix.
}
\usage{
rowMAD(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    \code{matrix}
  }
  \item{y}{
    ignored
  }
  \item{\dots}{
    Addition arguments to function \code{\link{mad}}.
  }
}

\value{
  A \code{numeric} vector of median absolute deviations.
}

\author{
R.Scharpf
}

\seealso{
  \code{\link{mad}}
}
\keyword{manip}

