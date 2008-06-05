\name{glass}
\alias{glass}


\docType{data}
\title{The Glass Data files for GAMLSS}
\description{

  glass: show the \code{strength} of glass fibres, measured at the National Physical Laboratory, England,
           see Smith  and Naylor (1987), (the unit of measurement were not given in the paper).
  }
\usage{
data(glass)
}
\format{
  Data frames each with the following variable.
  \describe{
    \item{\code{strength}}{a numeric vector showing the strength of glass fibres}
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{


Smith R. L. Naylor, J. C. (1987) A comparison of maximum likelihood and Bayesian estimators for the three-parameter
Weibull distributuion. \emph{Appl. Statist.} \bold{36}, 358-369
}

\examples{
data(glass)
with(glass, hist(strength))
}
\keyword{datasets}
