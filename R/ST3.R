# Monday, November 27, 2004 at 16:50 
# the first derivatives squares have been used here
#----------------------------------------------------------------------------------------
ST3 <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Skew t, type 3", substitute(mu.link), 
                           c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Skew t, type 3", substitute(sigma.link), 
                           c("inverse", "log", "identity"))
    vstats <- checklink(   "nu.link", "Skew t, type 3", substitute(nu.link),    
                           c("1/nu^2", "log", "identity"))
    tstats <- checklink(  "tau.link", "Skew t, type 3", substitute(tau.link),   
                            c("inverse", "log", "identity")) 
    structure(
          list(family = c("ST3", "Skew t, type 3"),
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
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
      su <- sqrt(1-u^2)
    dldm <- (1-nu/r+tau/2)*su*(1+u)-(1+nu/r+tau/2)*su*(1-u)
    dldm <- dldm/(sigma*sqrt(2*tau))
                    },
   d2ldm2 = function(){
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
      su <- sqrt(1-u^2)
    dldm <- (1-nu/r+tau/2)*su*(1+u)-(1+nu/r+tau/2)*su*(1-u)
    dldm <- dldm/(sigma*sqrt(2*tau))
   d2ldm2 <- -dldm*dldm
                      },
   dldd = function() {  
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
    dldd <- (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldd <- -(1/sigma) + dldd/(sigma*tau) 
      } ,
   d2ldd2 = function(){
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
    dldd <- (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldd <- -(1/sigma) + dldd/(sigma*tau) 
      d2ldd2 <- -dldd*dldd
                      },   
     dldv = function() { 
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
     dldv <- -digamma((1+nu/r)/tau) + digamma((1-nu/r)/tau)+log(1+u)-log(1-u)
     dldv <- 2*dldv/(r^3)
                        } ,
    d2ldv2 = function() { 
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
     dldv <- -digamma((1+nu/r)/tau) + digamma((1-nu/r)/tau) +log(1+u)-log(1-u)
     dldv <- 2*dldv/(r^3) 
      d2ldv2 <-  -dldv*dldv             
                        },
      dldt = function() {
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
       dldt <- 4*digamma(2/tau)-tau-4*log(2)
       dldt <- dldt + (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
       dldt <- dldt + 2*(1+nu*(r^2+tau)/r^3)*(log(1+u)-digamma((1+nu/r)/tau))
       dldt <- dldt + 2*(1-nu*(r^2+tau)/r^3)*(log(1-u)-digamma((1-nu/r)/tau))
       dldt <- -dldt/(2*tau*tau)
                        } ,
      d2ldt2 = function() { 
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
       dldt <- 4*digamma(2/tau)-tau-4*log(2)
       dldt <- dldt + (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
       dldt <- dldt + 2*(1+nu*(r^2+tau)/r^3)*(log(1+u)-digamma((1+nu/r)/tau))
       dldt <- dldt + 2*(1-nu*(r^2+tau)/r^3)*(log(1-u)-digamma((1-nu/r)/tau))
       dldt <- -dldt/(2*tau*tau)
      d2ldt2 <-   -dldt*dldt   
                            } ,
  d2ldmdd = function() {
         d2ldmdd <- -(dldm*dldd)
                       },
  d2ldmdv = function() {
         d2ldmdv <- -(dldm*dldv)
                       },
  d2ldmdt = function() {
         d2ldmdt <- -(dldm*dldt)
                       },
  d2ldddv = function() {
         d2ldddv <- -(dldd*dldv)
                       },
  d2ldddt = function() {
         d2ldddt <- -(dldd*dldt)  
                       },
  d2ldvdt = function() {
         d2ldvdt <- -(dldv*dldt)  
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
G.dev.incr <- -2*dST3(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(
            rqres(pfun="pST3", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)   
                           ) ,
    mu.initial = expression(mu <- (y+mean(y))/2) ,    #(y+mean(y))/2),# rep(mean(y),length(y)) 
 sigma.initial = expression(sigma <- rep(sd(y)/4, length(y))),
    nu.initial = expression(nu <- rep(0.03, length(y))), 
   tau.initial = expression(tau <-rep(3, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dST3 <- function(y, mu = 0, sigma = 1, nu = 0, tau = 1, log = FALSE)
 {
       if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
       if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
      z <- (y-mu)/sigma
     ve <- 2/tau
    lam <- 2*nu/(tau*sqrt(2*tau+nu*nu))
      a <- (ve+lam)/2
      b <- (ve-lam)/2
 loglik <- (a+0.5)*log(1+z/(sqrt(a+b+z*z)))+(b+0.5)*log(1-z/(sqrt(a+b+z*z)))-(a+b-1)*log(2)-0.5*log(a+b)-lbeta(a, b)-log(sigma)
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#----------------------------------------------------------------------------------------
pST3 <- function(q, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE)
 {  
      if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
      if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))           
      z <- (q-mu)/sigma
     ve <- 2/tau
    lam <- 2*nu/(tau*sqrt(2*tau+nu*nu))
      a <- (ve+lam)/2
      b <- (ve-lam)/2
  alpha <- (1+z/(sqrt(a+b+z*z)))/2
# we need the incomplete beta function ratio, i.e. cdf of Beta(a,b) at alpha
      p <- pbeta(alpha,a,b)
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#----------------------------------------------------------------------------------------
qST3 <-  function(p, mu=0, sigma=1, nu=0, tau=1, lower.tail = TRUE, log.p = FALSE)
 {   
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
# we need the inverse cdf of Beta(a,b) corresponding to probability p (i.e. our alpha)
    ve <- 2/tau
    lam <- 2*nu/(tau*sqrt(2*tau+nu*nu))
    a <- (ve+lam)/2
    b <- (ve-lam)/2
    Balpha <- qbeta(p,a,b) 
    zalpha <- (sqrt(a+b))*(2*Balpha-1)/(2*sqrt(Balpha*(1-Balpha)))
    q <- mu + sigma*zalpha  
    q
 }
#----------------------------------------------------------------------------------------
rST3 <- function(n, mu=0, sigma=1, nu=0, tau=1)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qST3(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
#----------------------------------------------------------------------------------------
