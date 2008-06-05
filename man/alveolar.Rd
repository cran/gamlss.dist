\name{alveolar}

\alias{alveolar}


\docType{data}
\title{The Alveolar Data files for GAMLSS}
\description{
 

  alveolar : alveolar-bronchiolar adenomas data used by Tamura and Young (1987) and also reproduce in Hand {\it et
  al.}   (1994), data set 256. The data are the number of mice out of certain number of mice (the binomial denominator) in
             23 independent  groups, having alveolar-bronchiolar adenomas.
}
\usage{


data(alveolar)

}
\format{
  Data frames each with the following variable.
  \describe{

  
  
    \item{\code{r}}{a numeric vector showing the number of mice out of
             n number of mice (the binomial denominator below) in 23 independent
             groups, having alveolar-bronchiolar adenomas.}
    \item{\code{n}}{a numeric vector showing the total number of mice}
    
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{


Hand \emph{et al.} (1994) \emph{A handbook of small data sets}. Chapman and Hall, London.


}

\examples{
data(alveolar)
with(alveolar, hist(r/n))

}
\keyword{datasets}
