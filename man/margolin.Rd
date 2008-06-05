\name{margolin}

\alias{margolin}


\docType{data}
\title{The Margolin Data files for GAMLSS}
\description{
 

  margolin: Margolin et al. (1981) present data from an Ames Salmonella assay,
  where y is the number of revertant colonies observed on a plate given a dose y of quinoline.
  The data were subsequently analysed by Breslow (1984), Lawless (1987) and Saha and Paul (2005).

 
  }
\usage{
data(parzen)
data(glass)
data(tensile)
data(margolin)
data(computer)
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
   
    \item{\code{y}}{a numeric vector showing the number of revertant colonies observed on a plate given a dose
                x of quinoline.}
    \item{\code{x}}{a numeric vector showing a a dose
                x of quinoline.}
   
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{
Breslow, N. (1984) Extra-Poisson variation in log-linear models. \emph{Applied Statistics}, \bold{33}, 38-44.


Hand \emph{et al.} (1994) \emph{A handbook of small data sets}. Chapman and Hall, London.

Lawless, J.F. (1987) Negative binomial and mixed Poisson regression. \emph{The Canadian Journal of Statistics},
\bold{15}, 209-225.

Margolin, B.H., Kaplan, N. and Zeiger, E. (1981) Statistical analysis of the Ames salmonella/microsome test.
  \emph{Proceedings of the National Academy of Science}, U.S.A., \bold{76}, 3779-3783.


Saha, K. and Paul, S. (2005) Bias-Corrected Maximum Likelihood Estimator of the Negative Binomial Dispersion
Parameter.
          \emph{Biometrics}, \bold{61}, 179-185

}

\examples{

data(margolin)
with(margolin, plot(y~x))

}
\keyword{datasets}
