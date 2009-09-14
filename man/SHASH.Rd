\name{SHASH}
\alias{SHASH}
\alias{dSHASH}
\alias{pSHASH}
\alias{qSHASH}
\alias{rSHASH}



\title{The Sinh-Arcsinh (SHASH) distribution for fitting a GAMLSS}
\description{
This function defines the Sinh-Arcsinh (SHASH) distribution, a four parameter distribution, 
for a \code{gamlss.family} object to be used for a  
GAMLSS fitting using the function \code{gamlss()}.  The functions \code{dSHASH}, 
\code{pSHASH}, \code{qSHASH} and \code{rSHASH} define the density, 
distribution function, quantile function and random
generation for the Sinh-Arcsinh (SHASH) distribution.   
}
\usage{
SHASH(mu.link = "identity", sigma.link = "log", nu.link = "log", 
      tau.link = "log")
dSHASH(x, mu = 0, sigma = 1, nu = 0.5, tau = 0.5, log = FALSE)
pSHASH(q, mu = 0, sigma = 1, nu = 0.5, tau = 0.5, lower.tail = TRUE, 
     log.p = FALSE)
qSHASH(p, mu = 0, sigma = 1, nu = 0.5, tau = 0.5, lower.tail = TRUE, 
     log.p = FALSE)
rSHASH(n, mu = 0, sigma = 1, nu = 0.5, tau = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{mu.link}{Defines the \code{mu.link}, with "identity" link as the default for the \code{mu} parameter. }
 \item{sigma.link}{  Defines the  \code{sigma.link}, with "log" link as the default for the \code{sigma} parameter. }
 \item{nu.link}{Defines the  \code{nu.link}, with "log" link as the default for the \code{nu} parameter.}
 \item{tau.link}{Defines the  \code{tau.link}, with "log" link as the default for the \code{tau} parameter.}
 \item{x,q}{vector of quantiles}
 \item{mu}{vector of location parameter values }
 \item{sigma}{vector of scale parameter values}
 \item{nu}{vector of skewness \code{nu} parameter values}
 \item{tau}{vector of kurtosis \code{tau} parameter values}
 \item{log, log.p}{ logical; if TRUE, probabilities p are given as log(p).}
 \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
          otherwise, P[X > x] }
 \item{p}{vector of probabilities. }
 \item{n}{ number of observations. If \code{length(n) > 1}, the length is
          taken to be the number required}      
}
\details{   
 The probability density function of the Sinh-Arcsinh distribution, (\code{SHASH}), Jones(2005), is defined as 
  
\deqn{f(y|\mu,\sigma\,\nu,\tau)= \frac{c}{\sqrt{2 \pi} \sigma (1+z^2)^{1/2}} e^{-r^2/2}}{f(y|mu,sigma,nu,tau)=c/(sqrt(2*pi)*sigma*(1+z^2)^(1/2)) exp(-(r^2)/2)}

where 

\deqn{r=\frac{1}{2} \left \{ \exp\left[ \tau \sinh^{-1}(z) \right] -\exp\left[ -\nu \sinh^{-1}(z)  \right] \right\}}{r=0.5*[exp(tau*sinh^(-1)(z))-exp(-nu*sinh^(-1)(z))]}

and

\deqn{c=\frac{1}{2} \left \{ \tau \exp\left[ \tau \sinh^{-1}(z) \right] + \nu \exp\left[ -\nu \sinh^{-1}(z)  \right] \right\}}{c=0.5*[rho*exp(tau*sinh^(-1)(z))+ nu * exp(-nu*sinh^(-1)(z))]}

and  \eqn{z=(y-\mu)/\sigma}{z=(y-mu)/sigma}
for \eqn{ -\infty < y < \infty }{0<y<0},
\eqn{\mu=(-\infty,+\infty)}{mu=(-Inf,+Inf)},
\eqn{\sigma>0}{sigma>0},
\eqn{\nu>0}{nu>0} and
\eqn{\tau>0}{tau>0}. 



The parameters \eqn{\mu}{mu} and \eqn{\sigma}{sigma} are the location  and scale of the distribution.
The parameter \eqn{\nu}{nu} determines the left hand tail of the distribution with \eqn{\nu>1}{nu>1} indicating a lighter tail  than the normal 
and 
\eqn{\nu<1} heavier tail than the normal. The parameter \eqn{\tau}{tau}  determines the right hand tail  of the distribution in the same way.

}
\value{
\code{SHASH()}  returns a \code{gamlss.family} object which can be used to fit the SHASH distribution in the \code{gamlss()} function. 
\code{dSHASH()} gives the density, \code{pSHASH()} gives the distribution
     function, \code{qSHASH()} gives the quantile function, and \code{rSHASH()}
     generates random deviates. 
}
\references{

Jones, M. C. (2006) p 546-547 in the discussion of Rigby, R. A. and  Stasinopoulos D. M. (2005) 
\emph{Appl. Statist.}, \bold{54}, part 3.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R. A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).  

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

}
\author{Bob Rigby \email{r.rigby@londonmet.ac.uk} and Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }


\section{Warning}{The qSHASH and rSHASH are slow since they are relying on golden section for finding the quantiles} 

\seealso{  \code{\link{gamlss.family}}, \code{\link{JSU}}, \code{\link{BCT}}}
\examples{
SHASH()   # 
plot(function(x)dSHASH(x, mu=0,sigma=1, nu=1, tau=2), -5, 5, 
 main = "The SHASH  density mu=0,sigma=1,nu=1, tau=2")
plot(function(x) pSHASH(x, mu=0,sigma=1,nu=1, tau=2), -5, 5, 
 main = "The BCPE  cdf mu=0, sigma=1, nu=1, tau=2")
dat<-rSHASH(100,mu=10,sigma=1,nu=1,tau=1.5)
hist(dat)
# library(gamlss)
# gamlss(dat~1,family=SHASH, control=gamlss.control(n.cyc=30))
}
\keyword{distribution}
\keyword{ regression}
