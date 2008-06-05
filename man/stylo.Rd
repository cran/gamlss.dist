\name{stylo}
\alias{stylo}

\docType{data}
\title{The Stylometric Data files for GAMLSS}
\description{
 stylo : the data were given by Dr Mario Corina-Borja, see Chappas and Corina-Borja (2006), and has the number of a
 word appearing in a text.
  }
\usage{
data(stylo)
}
\format{
  Data frames each with the following variable.
  \describe{

    \item{\code{word}}{a numeric vector showing the number  a word appearing in a text}
    \item{\code{freq}}{a numeric vector showing the frequency of the number  a word appearing in a text}
  }
}
\details{Data sets usefull for the GAMLSS booklet}
\references{
Chappas C. and Corina-Borja M. A Stylometric analysis of newspapers periodical and news scriprs, 
\emph{Journal of Quantitative Linguistics}, 13, 285-312}

\examples{
data(stylo)
plot(freq~word, type="h", data=stylo)
}
\keyword{datasets}
