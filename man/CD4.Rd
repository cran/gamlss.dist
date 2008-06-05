\name{CD4}

\alias{CD4}


\docType{data}
\title{The CD4 Count Data files for GAMLSS}
\description{

 CD4: The data were given by Wade and Ades (1994) and refer to cd4 counts from uninfected children
   born to HIV-1 mothers and the age of the child.


  }
\usage{

data(CD4)
}
\format{
  Data frames each with the following variable.
  \describe{




   
    \item{cd4}{a numeric vector showing the CD4 counts}
    \item{age}{the age of the child}
  
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{
Wade, A. M. and Ader, A. E. (1994) Age-related reference ranges : Significance tests for models and confidence
intervals for centiles. \emph{Statistics in Medicine}, \bold{13}, pages 2359-2367.
}

\examples{



data(CD4)
with(CD4,plot(cd4~age))

}
\keyword{datasets}
