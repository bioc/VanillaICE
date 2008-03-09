\name{makeTable}
\alias{makeTable}
\title{Makes a table of summary statistics for the breakpoints}
\description{
  Calculate summary statistics of the breakpoints (mean, median, sd) for
  each chromosome.  The average, median, and sd are also computed across chromosomes.
}
\usage{
makeTable(object, state, by, unit = 1000, digits = 3)
}
\arguments{
  \item{object}{ object of class \code{HmmPredict}}
  \item{state}{a character string for the hidden state to make the
    summary table of}
  \item{by}{not useful at this point.  Eventually would be nice to split
  by sample as well as by chromosome}
  \item{unit}{numeric.  For kb, set unit = 1000.  For Mb, set unit =
    1e6, etc.}
  \item{digits}{for rounding}
}
\value{
  matrix
}
\author{R. Scharpf}
\seealso{\code{\link{summary}}}
\examples{
data(fit)
makeTable(fit, state="D")
}
\keyword{manip}

