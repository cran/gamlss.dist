\name{computer}

\alias{computer}

\docType{data}
\title{The Computer Failure Data files for GAMLSS}
\description{
 

  computing: The data relate to DEC-20 computers which operated at the Open University in the 1980. They give the
  number of
  computers that broke down in each of the 128 consecutive weeks of operation, starting in late 1983, see Hand {\it et
  al.} (1994) page 109 data set 141.

  }
\usage{

data(computer)

}
\format{
  Data frames each with the following variable.
  \describe{
    \item{\code{failure}}{a numeric vector showing the number of times computers failed}
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{

Hand \emph{et al.} (1994) \emph{A handbook of small data sets}. Chapman and Hall, London.


}

\examples{
data(computer)
with(computer, plot(table(failure)))

}
\keyword{datasets}
