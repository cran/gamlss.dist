# Amended 20/10/2007, see change to sigma link, dEGB2 and dldd to allow for negative sigma
# starting values for mu and sigma  (and nu and tau) may need to be changed 
EGB2 <- function (mu.link="identity", sigma.link="identity", nu.link ="log", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Exponential generalized beta 2 (i.e. of the second kind)", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Exponential generalized beta 2 (i.e. of the second kind)", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Exponential generalized beta 2 (i.e. of the second kind)", substitute(nu.link),    
                           c("1/nu^2", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Exponential generalized beta 2 (i.e. of the second kind)", substitute(tau.link),   
                           c("1/tau^2", "log", "identity", "own")) 
    structure(
          list(family = c("EGB2", "Exponential generalized beta 2 (i.e. of the second kind)"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
             tau.link = as.character(substitute(tau.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
           tau.linkfun = tstats$linkfun,  
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
           tau.linkinv = tstats$linkinv, 
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
               tau.dr = tstats$mu.eta, 
    dldm = function() { 
      z <- (y-mu)/sigma
   dldm <- ((exp(z)*(tau+nu))/(sigma*(1+exp(z))))-(nu/sigma)
                      },
   d2ldm2 = function(){
      d2ldm2 <- -(nu*tau)/((nu+tau+1)*sigma^2)
      d2ldm2
                      },     
   dldd = function() {  
      z <- (y-mu)/sigma
      dldd <- z*((exp(z)*(tau+nu))/(sigma*(1+exp(z))))-(nu*z/sigma)-1/sigma
      } ,
   d2ldd2 = function(){
      d2a <- nu*(digamma(nu+1)-digamma(nu))
      d2b <- ((nu*tau)/(nu+tau+1))*(trigamma(nu+1)+trigamma(tau+1)+(digamma(nu)-digamma(tau))^2)
      d2ldd2 <- -(d2a+d2b)/(sigma^2)
      d2ldd2
                      },   
     dldv = function() { 
      z <- (y-mu)/sigma
      dldv <- z-log(1+exp(z))-digamma(nu)+digamma(nu+tau)
                        } ,
    d2ldv2 = function() { 
      d2ldv2 <- trigamma(nu+tau) - trigamma(nu)
      d2ldv2
                        },
      dldt = function() {
      z <- (y-mu)/sigma  
      dldt <- -log(1+exp(z))-digamma(tau)+digamma(nu+tau)
                        } ,
      d2ldt2 = function() { 
      d2ldt2 <- trigamma(nu+tau) - trigamma(tau)
      d2ldt2
                            } ,
  d2ldmdd = function() {
  d2ldmdd <- -((nu*tau)/((nu+tau+1)*sigma^2))*(digamma(nu+1)-digamma(tau+1))
  d2ldmdd 
                       },
  d2ldmdv = function() {
  d2ldmdv <- -tau/(sigma*(nu+tau))
  d2ldmdv
                       },
  d2ldmdt = function() {
  d2ldmdt <- nu/(sigma*(nu+tau))
  d2ldmdt
                       },
  d2ldddv = function() {
  d2ldddv <- nu/(sigma*(nu+tau))*(digamma(nu+1)-digamma(tau))-(1/sigma)*(digamma(nu)-digamma(tau))
  d2ldddv
                       },
  d2ldddt = function() {
  d2ldddt <- (nu/(sigma*(nu+tau)))*(digamma(nu+1)-digamma(tau))
  d2ldddt
                       },
  d2ldvdt = function() {
  d2ldvdt <- trigamma(nu+tau)
  d2ldvdt
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
G.dev.incr <- -2*dEGB2(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               {                                                               
               cdf <- pEGB2(y,mu,sigma,nu,tau)           
             rqres <- qnorm(cdf)
                }) ,
    mu.initial = expression(mu <- (y+mean(y))/2), 
    sigma.initial = expression(sigma <- rep(0.1, length(y))), 
    nu.initial = expression(nu <- rep(1, length(y))), 
   tau.initial = expression(tau <-rep(1, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma) TRUE,
      nu.valid = function(nu) all(nu > 0), 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dEGB2 <- function(y, mu = 0, sigma = 1, nu = 1, tau = .5, log = FALSE)
 {
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
      z <- (y-mu)/sigma
 loglik <- nu*z-log(abs(sigma))-lgamma(nu)-lgamma(tau)+lgamma(nu+tau)-(nu+tau)*log(1+exp(z))
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pEGB2 <- function(q, mu = 0, sigma = 1, nu = 1, tau = .5, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))           
      w <- (tau/nu)*exp((q-mu)/sigma)  
      p <- pf(w,2*nu,2*tau)
    if (length(sigma)>1)  p <- ifelse(sigma<0,1-p,p)
    else p <- if (sigma<0) 1-p else p
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#-----------------------------------------------------------------  

qEGB2 <-  function(p, mu=0, sigma=1, nu=0, tau=.5, lower.tail = TRUE, log.p = FALSE)
 {   
    if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (length(sigma)>1)  p <- ifelse(sigma<0,1-p,p)
    else p <- if (sigma<0) 1-p else p
    w <- qf(p,2*nu,2*tau)   
    q <- mu+sigma*log((nu/tau)*w)  
    q
 }
#-----------------------------------------------------------------  
rEGB2 <- function(n, mu=0, sigma=1, nu=0, tau=.5)
  {
    if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qEGB2(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
#-----------------------------------------------------------------  
