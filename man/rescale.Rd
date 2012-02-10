\name{rescale}
\alias{rescale}
\title{Rescale a numeric vector}
\description{
  Rescale a numeric vector
}
\usage{
rescale(x, l, u)
}

\arguments{
  \item{x}{
    a \code{numeric} vector
}
  \item{l}{
    \code{numeric}: lower limit of rescaled \code{x}.
}
  \item{u}{
    \code{numeric}: upper limit of rescaled \code{x}.
}
}
\details{
  Not intended to be called directly, but used in packages that depend
  on \pkg{VanillaICE}
}
\value{
  \code{numeric} vector the same length as \code{x} with range [l, u].
}

\author{
R. Scharpf
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
