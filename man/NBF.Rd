\name{NBF}
\alias{NBF}
\alias{dNBF}
\alias{pNBF}
\alias{qNBF}
\alias{rNBF}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Negative Binomial Family distribution for fitting a GAMLSS}
\description{
The \code{NBF()} function defines the Negative Binomial family distribution, a three parameter distribution, for a \code{gamlss.family} object to be used  in GAMLSS fitting using the function \code{gamlss()}.    
The functions \code{dNBF}, \code{pNBF}, \code{qNBF} and \code{rNBF} define the density, distribution function, quantile function and random
generation for the Negative Binomial family, \code{NBF()}, distribution
}
\usage{
NBF(mu.link = "log", sigma.link = "log", nu.link = "identity")

dNBF(x, mu = 0, sigma = 1, nu = 2, log = FALSE)

pNBF(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

qNBF(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

rNBF(n, mu = 0, sigma = 1, nu = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu.link}{The link function for \code{mu}}
  \item{sigma.link}{The link function for \code{sigma}}
  \item{nu.link}{The link function for \code{nu}}
  \item{x}{vector of (non-negative integer)}
  \item{mu}{vector of positive means}
  \item{sigma}{vector of positive despersion parameter}
  \item{nu}{vector of power parameter}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p)}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]}
  \item{p}{vector of probabilities}
  \item{q}{vector of quantiles}
  \item{n}{number of random values to return}
}

\details{
 The definition for Negative Binomial Family distribution , \code{NBF}, is similar to the Negative Binomial type I. The probability function of the \code{NBF} can be obtained by replacing \eqn{\sigma} with  \eqn{\sigma \mu^{\nu-2}} where \eqn{\nu} is a power parameter.
The distrobution has mean \eqn{\mu} and variance \eqn{\mu+\sigma \mu^{vu}.} 
}
\value{
  returns a \code{gamlss.family} object which can be used to fit a Negative Binomial Family distribution in the \code{gamlss()} function. 
}
\references{
Anscombe, F. J. (1950) Sampling theory of the negative bimomial and logarithmic distributiona, \emph{Biometrika}, \bold{37}, 358-382.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).  

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Bob Rigby and Mikis Stasinopoulos}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{NBI}}, \code{\link{NBII}}}
\examples{
NBF()   # gives information about the default links for the Negative Binomial Family 
# plotting the distribution
plot(function(y) dNBF(y, mu = 10, sigma = 0.5, nu=2 ), from=0, to=40, n=40+1, type="h")
# creating random variables and plot them 
tN <- table(Ni <- rNBF(1000, mu=5, sigma=0.5, nu=2))
r <- barplot(tN, col='lightblue')
\dontrun{
library(gamlss)
data(species)
species <- transform(species, x=log(lake))
m6 <- gamlss(fish~poly(x,2), sigma.fo=~1, data=species, family=NBF, n.cyc=200)
fitted(m6, "nu")[1]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{regression}
