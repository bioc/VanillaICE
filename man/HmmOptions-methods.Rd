\name{HmmOptions-methods}
\alias{HmmOptions-methods}
\alias{copyNumber.location}
\alias{copyNumber.location<-}
\alias{probHomCall}
\alias{probHomCall<-}
\alias{snpset}
\alias{snpset<-}
\alias{term5}
\alias{term5<-}
\title{Accessors and replacement methods for HmmOptions}
\description{Accessors and replacement methods for HmmOptions objects}
\usage{

	copyNumber.location(object)
	copyNumber.location(object) <- value	
	probHomCall(object)
	probHomCall(object) <- value
	snpset(object)
	snpset(object) <- value
	term5(object)
	term5(object) <- value
}

\arguments{
	\item{object}{\code{HmmOptions} instance}
	\item{value}{see the signature of the corresponding methods}
      }

\value{

  copyNumber.location, probHomCall and term5 each return numeric vectors.


  snpset returns an instance of \code{SnpLevelSet}
  
}
\seealso{\code{\link{HmmOptions-class}}}
\keyword{methods}
       

