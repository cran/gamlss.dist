\name{parzen}
\alias{parzen}
\alias{glass}
\alias{tensile}
\alias{margolin}
\alias{computer}
\alias{lice}
\alias{alveolar}
\alias{species}
\alias{CD4}
\alias{LGAclaims}

\docType{data}
\title{Data files for GAMLSS}
\description{
  These are several small data files usefull for gamlss fits.
  
  Parzen: Parzen (1979) and also contained in Hand {\it et al.} (1994), data set 278. The data
         give the annual \code{snowfall} in Buffalo, NY (inches) for the 63 years, from 1910 to 1972 inclusive.
  
  glass: show the \code{strength} of glass fibres, measured at the National Physical Laboratory, England, 
           see Smith  and Naylor (1987), (the unit of measurement were not given in the paper).  
  
  tensile: These data  come from Quesenberry and Hales (1980)  and were also reproduced in Hand {\it et al.} (1994),
    data set 180, page 140. They contain measurements of tensile strength of polyester fibres and the authors 
    were trying to check if they were consistent with the lognormal distribution. 
    According to Hand {\it et al.} (1994) "these data follow from a preliminary transformation.
    If the lognormal hypothesis is correct, these data should have been uniformly distributed".
  
  margolin: Margolin et al. (1981) present data from an Ames Salmonella assay, 
  where y is the number of revertant colonies observed on a plate given a dose y of quinoline.    
  The data were subsequently analysed by Breslow (1984), Lawless (1987) and Saha and Paul (2005).
 
  computing: The data relate to DEC-20 computers which operated 
at the Open University in the 1980. They give the number of
computers that broke down in each of the 128 consecutive weeks of 
operation, starting in late 1983, see Hand {\it et al.} (1994) page 109 data set 141.  
    
  lice : The data come from Williams (1944) and they are lice per head of Hindu male prisoners in 
              Cannamore, South India, 1937-1939.  
              
  alveolar : alveolar-bronchiolar adenomas data used by Tamura and Young (1987) and also
reproduce in Hand {\it et al.} (1994), data set 256. The data are the number of mice out of 
certain number of mice (the binomial denominator) in 23 independent
groups, having alveolar-bronchiolar adenomas.

 species: The number of different fish species (\code{y=fish}) was recorded
for 70 lakes of the world together with explanatory variable
\code{x=log(lake)} area.  The data are given and analyzed by Stein and Juritz (1988). 

 CD4: The data were given by Wade and Ades (1994) and refer to cd4 counts from uninfected children 
   born to HIV-1 mothers and the age of the child.
 
 LGAclaims: the data were given by Gillian Heller and can be found in de Jong and Heller (2007). 
This data set records the number of third party claims, \code{Claims}, in a twelve month
period between 1984-1986 in each of 176 geographical areas (local government areas) in New South Wales, 
 Australia. Areas are grouped into thirteen statistical divisions (\code{SD}). Other 
 recorded variables are the number of accidents, 
 \code{Accidents},  the number of people killed or injured and population with all variables classified
 according to area. 

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
}
\format{
  Data frames each with the following variable.
  \describe{
    \item{\code{snowfall}}{the annual snowfall in Buffalo, NY (inches) for the 63 years, from 1910 to 1972 inclusive, 63 observations}
    \item{\code{strength}}{a numeric vector showing the strength of glass fibres}
    \item{\code{str}}{a numeric vector showing the tensile strength}
    \item{\code{y}}{a numeric vector showing the number of revertant colonies observed on a plate given a dose 
                x of quinoline.}
    \item{\code{x}}{a numeric vector showing a a dose 
                x of quinoline.}
    \item{\code{failure}}{a numeric vector showing the number of times computers failed}
    \item{\code{head}}{a numeric vector showing the number lice per head of Hindu male prisoners in 
              Cannamore, South India, 1937-1939.}
    \item{\code{freq}}{a numeric vector showing the frequency of lice per head}
    \item{\code{r}}{a numeric vector showing the number of mice out of 
             n number of mice (the binomial denominator below) in 23 independent
             groups, having alveolar-bronchiolar adenomas.}
    \item{\code{n}}{a numeric vector showing the total number of mice}
    \item{\code{fish}}{a numeric vector showing the number of different species in 70 lakes in the word}
    \item{\code{lake}}{a numeric vector showing the lake area}
    \item{cd4}{a numeric vector showing the CD4 counts}
    \item{age}{the age of the child}
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
Breslow, N. (1984) Extra-Poisson variation in log-linear models. \emph{Applied Statistics}, \bold{33}, 38-44.

de Jong, P. and Heller G. (2007)  \emph{Generalized Linear Models for Insurance Data }, Cambridge University Press

Hand \emph{et al.} (1994) \emph{A handbook of small data sets}. Chapman and Hall, London. 

Lawless, J.F. (1987) Negative binomial and mixed Poisson regression. \emph{The Canadian Journal of Statistics}, \bold{15}, 209-225.

Margolin, B.H., Kaplan, N. and Zeiger, E. (1981) Statistical analysis of the Ames salmonella/microsome test. 
  \emph{Proceedings of the National Academy of Science}, U.S.A., \bold{76}, 3779-3783. 

Quesenberry, C. and Hales, C. (1980). Concentration bands for uniformily plots.
\emph{Journal of Statistical Computation and Simulation}, \bold{11}, 41:53.

Parzen E. (1984) Nonparamemetric statistical daya modelling. \emph{JASA}, \bold{74}, 105-131.

Saha, K. and Paul, S. (2005) Bias-Corrected Maximum Likelihood Estimator of the Negative Binomial Dispersion Parameter. 
          \emph{Biometrics}, \bold{61}, 179-185

Smith R. L. Naylor, J. C. (1987) A comparison of maximum likelihood and Bayesian estimators for the three-parameter
Weibull distributuion. \emph{Appl. Statist.} \bold{36}, 358-369

Stein, G. Z. and Juritz, J. M. (1988). Linear models with an inverse Gaussian-Poisson
error distribution. \emph{Communications in Statistics- Theory and Methods}, \bold{17}, 557-571.

Wade, A. M. and Ader, A. E. (1994) 
Age-related reference ranges : Significance tests for models and confidence 
intervals for centiles. \emph{Statistics in Medicine}, \bold{13}, pages 2359-2367.

}
\examples{
data(parzen)
with(parzen, hist(snowfall))
data(glass)
with(glass, hist(strength))
data(tensile)
with(tensile,hist(str))
data(margolin)
with(margolin, plot(y~x))
data(computer)
with(computer, plot(table(failure)))
data(lice)
with(lice, plot(freq~head, type="h"))
data(alveolar)
with(alveolar, hist(r/n))
data(species)
with(species, plot(fish~log(lake)))
data(CD4)
with(CD4,plot(cd4~age))
data(LGAclaims)
with(LGAclaims, plot(data.frame(Claims, Pop_density, KI, Accidents, Population)))
}
\keyword{datasets}
