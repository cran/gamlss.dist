\name{lice}

\alias{lice}


\docType{data}
\title{Data files for GAMLSS}
\description{
 
  lice : The data come from Williams (1944) (also used by Stein and Juritz (1988).) and they are lice per head of Hindu male prisoners in
              Cannamore, South India, 1937-1939.

 
  }
\usage{

data(lice)
data(alveolar)
data(species)
data(CD4)
data(LGAclaims)
data(stylo)
}
\format{
  Data frames each with the following variable.
  \describe{

    \item{\code{head}}{a numeric vector showing the number lice per head of Hindu male prisoners in
              Cannamore, South India, 1937-1939.}
    \item{\code{freq}}{a numeric vector showing the frequency of lice per head}
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{

Stein, G. Z. and Juritz, J. M. (1988). Linear models with an inverse Gaussian-Poisson
error distribution. \emph{Communications in Statistics- Theory and Methods}, \bold{17}, 557-571.


}

\examples{

data(lice)
}
\keyword{datasets}
