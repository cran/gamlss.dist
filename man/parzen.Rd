\name{parzen}
\alias{parzen}


\docType{data}
\title{The Parzen Data File for GAMLSS}
\description{
    
  Parzen: Parzen (1979) and also contained in Hand \emph{et al.} (1994), data set 278. The data
         give the annual \code{snowfall} in Buffalo, NY (inches) for the 63 years, from 1910 to 1972 inclusive.
  
}
\usage{
data(parzen)
}
\format{
  Data frames each with the following variable.
  \describe{
    \item{\code{snowfall}}{the annual snowfall in Buffalo, NY (inches) for the 63 years, from 1910 to 1972 inclusive, 63 observations} 
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{

Hand \emph{et al.} (1994) \emph{A handbook of small data sets}. Chapman and Hall, London. 

Parzen E. (1984) Nonparamemetric statistical daya modelling. \emph{JASA}, \bold{74}, 105-131.
}

\examples{
data(parzen)
with(parzen, hist(snowfall))
}
\keyword{datasets}
