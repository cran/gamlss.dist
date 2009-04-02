\name{tensile}

\alias{tensile}

\docType{data}
\title{The Tensile Data files for GAMLSS}
\description{
 

  tensile: These data  come from Quesenberry and Hales (1980)  and were also reproduced in Hand \emph{et al.} (1994),
    data set 180, page 140. They contain measurements of tensile strength of polyester fibres and the authors
    were trying to check if they were consistent with the lognormal distribution.
    According to Hand \emph{et al.} (1994) "these data follow from a preliminary transformation.
    If the lognormal hypothesis is correct, these data should have been uniformly distributed".

 
  }
\usage{

data(tensile)

}
\format{
  Data frames each with the following variable.
  \describe{
    
    \item{\code{str}}{a numeric vector showing the tensile strength}
    
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{

Hand \emph{et al.} (1994) \emph{A handbook of small data sets}. Chapman and Hall, London.


Quesenberry, C. and Hales, C. (1980). Concentration bands for uniformily plots.
\emph{Journal of Statistical Computation and Simulation}, \bold{11}, 41:53.

}

\examples{

data(tensile)
with(tensile,hist(str))

}
\keyword{datasets}
