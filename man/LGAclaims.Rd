\name{LGAclaims}
\alias{LGAclaims}

\docType{data}
\title{The LGA Claims Data files for GAMLSS}
\description{
  These are several small data files usefull for gamlss fits.


 LGAclaims: the data were given by Gillian Heller and can be found in de Jong and Heller (2007).
       This data set records the number of third party claims, \code{Claims}, in a twelve month
      period between 1984-1986 in each of 176 geographical areas (local government areas) in New South Wales,
      Australia. Areas are grouped into thirteen statistical divisions (\code{SD}). Other
         recorded variables are the number of accidents, \code{Accidents},  the number of people killed or
         injured and population with all variables classified according to area.

  }
\usage{



data(LGAclaims)
}
\format{
  Data frames each with the following variable.
  \describe{
    \item{Claims}{the number of third party claims}
    \item{LGA}{Local government areas in New South Wales}
    \item{SD}{statistical divisions}
    \item{Pop$\_$density}{population density}
    \item{KI}{the number of people killed or injured}
    \item{Accidents}{the number of accidents}
    \item{Population}{population size}
    \item{L$\_$KI}{log of KI}
    \item{L$\_$Accidents}{the log of the number of accidents}
    \item{L$\_$Population}{log Population}
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{


de Jong, P. and Heller G. (2007)  \emph{Generalized Linear Models for Insurance Data }, Cambridge University Press


}

\examples{



data(LGAclaims)
with(LGAclaims, plot(data.frame(Claims, Pop_density, KI, Accidents, Population)))

}
\keyword{datasets}
