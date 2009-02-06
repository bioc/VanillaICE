\name{affy6}
\alias{affy6}
\docType{data}
\title{CRLMM confidence scores for affy 6.0}
\description{
  Confidence scores for affy 6.0 were generated using version
  1.0.34 of the crlmm package.
}
\usage{data(affy6)}
\format{
  The format is:
 int [1:179154, 1:2] 3 3 2 3 3 3 3 2 3 3 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "group" "phat"
}
\details{
  The phat confidence scores are stored as an integer to save memory.
  The confidence scores can be converted to a probability by

  prob <- 1-exp(-affy6[, "phat"]/1000)
}
\examples{
data(affy6)
}
\keyword{datasets}
