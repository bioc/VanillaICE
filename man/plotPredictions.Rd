\name{plotPredictions}
\alias{plotPredictions}
\title{Plots the predictions from the HMM}
\description{
Plots the predictions from the HMM.  This function is usually not called directly.
}
\usage{
plotPredictions(predictions, op, breakpoints, position, chromosome, states, ...)
}
\arguments{
  \item{predictions}{vector of predictions}
  \item{op}{an instance of a class for graphical parameters}
  \item{breakpoints}{matrix of breakpoints}
  \item{position}{vector of physical positions}
  \item{chromosome}{ vector of chromosomes}
  \item{states}{the hidden states}
  \item{\dots}{ not implemented}
}
\value{
  Nothing
}
\author{R. Scharpf}
\examples{
##See the VanillaICE vignette
}
\keyword{aplot}

