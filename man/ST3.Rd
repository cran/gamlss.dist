\name{ST3}
\alias{ST3}
\alias{dST3}
\alias{pST3}
\alias{qST3}
\alias{rST3}


\title{The skew t distribution, type 3 }
\description{
The Skew t type 3 distribution Jones and Faddy (2003). 
The functions \code{dST3}, \code{pST3}, \code{qST3} and \code{rST3} define the density, distribution function,
 quantile function and random generation for the skew t distribution type 3. 
      
}
\usage{
ST3(mu.link = "identity", sigma.link = "log", nu.link = "identity", tau.link="log")
dST3(y, mu = 0, sigma = 1, nu = 0, tau=1, log = FALSE)
pST3(q, mu = 0, sigma = 1, nu = 0, tau=1, lower.tail = TRUE, log.p = FALSE)
qST3(p, mu = 0, sigma = 1, nu = 0, tau=1, lower.tail = TRUE, log.p = FALSE)
rST3(n, mu = 0, sigma = 1, nu = 0, tau=1)
}

\arguments{
 \item{mu.link}{ Defines the \code{mu.link}, with "identity" link as the default for the \code{mu} parameter. 
                  Other links are "\eqn{1/mu^2}" and "log"}
  \item{sigma.link}{ Defines the  \code{sigma.link}, with "log" link as the default for the \code{sigma} parameter. 
                      Other links are "inverse" and "identity"}
  \item{nu.link}{Defines the  \code{nu.link}, with "identity" link as the default for the \code{nu} parameter. 
                  Other links are "\eqn{1/mu^2}" and "log"}
   \item{tau.link}{Defines the  \code{nu.link}, with "log" link as the default for the \code{nu} parameter. 
                   Other links are "inverse", "identity" }                
   \item{y,q}{vector of quantiles}
  \item{mu}{vector of \code{mu} parameter values }
  \item{sigma}{vector of scale parameter values}
  \item{nu}{vector of \code{nu} parameter values}
   \item{tau}{vector of \code{tau} parameter values}
  \item{log, log.p}{ logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
          otherwise, P[X > x] }
  \item{p}{vector of probabilities. }
  \item{n}{ number of observations. If \code{length(n) > 1}, the length is
          taken to be the number required}
  \item{\dots}{for extra arguments}
}
\details{

The probability density function of the skew t distribution type 3, (\code{ST3}), is defined as 
  
\deqn{f(y|\mu,\sigma,\nu, \tau)=\frac{1}{c} \left[ 1+ \frac{z}{(a+b +z^2)^{1/2}} \right]^{a+1/2} \left[ 1- \frac{z}{(a+b+z^2)^{1/2}}\right]^{b+1/2}}{f(y|mu,sigma,nu,tau)=(1/c)*(1+(z/(a+b+z^2)^0.5))^(a+0.5)*(1-(a+b+z^2)^0.5)^(b+0.5)}

where \eqn{c=z^{a +b-1} (a+b)^{1/2} B(a,b)}{c=z^(a+b-1) * (a+b)^0.5 *B(a,b)}, and 
\eqn{B(a,b)=\Gamma(a)\Gamma(b)/ \Gamma(a+b)}{Gamma(a)*Gamma(b)/Gamma(a+b)} and 
\eqn{z=(y-\mu)/\sigma}{(y-mu)/sigma} and 
\eqn{\nu=(a-b)/\left[ab(a+b) \right]^{1/2}}{nu=(a-b)/(a*b*(a+b))^0.5}
and
\eqn{\tau=2/(a+b)}{tau=2/(a+b)} 
for \eqn{-\infty<y>\infty}{-Inf<y<Inf}, \eqn{-\infty<\mu>\infty}{-Inf<mu<Inf}, \eqn{\sigma>0}, 
 \eqn{-\infty<\nu>\infty}{-Inf<nu<Inf} and \eqn{\tau>0}{tau>0}.   
}
\value{
\code{ST3()}  returns a \code{gamlss.family} object which can be used to fit the skew t type 3 distribution in the \code{gamlss()} function. 
  \code{dST3()} gives the density, \code{pST3()} gives the distribution  function, 
  \code{qST3()} gives the quantile function, and \code{rST3()}
  generates random deviates. 
  
}
\references{ 

Jones, M.C. and Faddy, M. J. (2003) A skew extension of the t distribution, with applications. 
\emph{Journal of the Royal Statistical Society}, Series B, \bold{65}, pp 159-174.  

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R. A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.londonmet.ac.uk/gamlss/}).

 }
\author{Bob Rigby and Mikis Stasinopoulos  }
\note{ The mean of the ex-Gaussian is \eqn{\mu+\nu}{mu+nu} and the variance is \eqn{\sigma^2+\nu^2}{sigma^2+nu^2}.
}
\seealso{ \code{\link[gamlss]{gamlss}}, \code{\link[gamlss]{gamlss.family}}, \code{\link[gamlss]{BCCG}}, \code{\link[gamlss]{GA}}, \code{\link[gamlss]{IG}}  \code{\link[gamlss]{LNO}}}
\examples{
%ST3()   # 
y<- rST3(100, mu=10, nu=1, sigma=3)
hist(y)
m1<-gamlss(y~1, family=ST3) 
plot(m1)
curve(dST3(y=x, mu=300 ,sigma=35,nu=100), 100, 600, 
 main = "The ex-GAUS  density mu=300 ,sigma=35,nu=100")
plot(function(x) pST3(x, mu=300,sigma=35,nu=100), 100, 600, 
 main = "The ex-GAUS  cdf mu=300, sigma=35, nu=100")
}
\keyword{distribution}% 
\keyword{regression}% 
