\name{species}


\alias{species}


\docType{data}
\title{The Fish Species Data files for GAMLSS}
\description{
 
 species: The number of different fish species (\code{y=fish}) was recorded for 70 lakes of the world together with 
 explanatory variable \code{x=log(lake)} area.  The data are given and analyzed by Stein and Juritz (1988).

   }
\usage{


data(species)

}
\format{
  Data frames each with the following variable.
  \describe{




    \item{\code{fish}}{a numeric vector showing the number of different species in 70 lakes in the word}
    \item{\code{lake}}{a numeric vector showing the lake area}
    
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{


Stein, G. Z. and Juritz, J. M. (1988). Linear models with an inverse Gaussian-Poisson
error distribution. \emph{Communications in Statistics- Theory and Methods}, \bold{17}, 557-571.

}

\examples{


data(species)
with(species, plot(fish~log(lake)))

}
\keyword{datasets}
