# the ex-Gaussian distribution
# Mikis Stasinopoulos and Bob Rigby (suggested by Jonathan Williams) 
# 21-3-06
exGAUS <- function (mu.link="identity", sigma.link="log", nu.link ="log") 
{
    mstats <- checklink("mu.link", "ex-Gaussian", substitute(mu.link), 
                         c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "ex-Gaussian", substitute(sigma.link), #
                         c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "ex-Gaussian",substitute(nu.link), 
                        c("logshifted", "log", "identity"), par.link = c(1))
    
    structure(
          list(family = c("exGAUS", "ex-Gaussian"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
                 dldm = function() {
        z <- y-mu-((sigma^2)/nu)
     dldm <- 1/nu-(dnorm(z/sigma)/(sigma*pnorm(z/sigma)))
     dldm
                                    },
               d2ldm2 = function() {
        z <- y-mu-((sigma^2)/nu)
     dldm <- 1/nu-(dnorm(z/sigma)/(sigma*pnorm(z/sigma)))
   d2ldm2 <- -dldm*dldm
                                    },
                 dldd = function() {
      z <- y-mu-((sigma^2)/nu) 
   pphi <- (dnorm(z/sigma)/(pnorm(z/sigma)))
   dldd <- (sigma/(nu^2))-((z/sigma^2)+(2/nu))*pphi
   dldd
                                    },
               d2ldd2 = function() {
      z <- y-mu-((sigma^2)/nu) 
   pphi <- (dnorm(z/sigma)/(pnorm(z/sigma)))
   dldd <- (sigma/(nu^2))-((z/sigma^2)+(2/nu))*pphi
 d2ldd2 <- -dldd*dldd
                                    },
                 dldv = function() {
       z <- y-mu-((sigma^2)/nu)
    pphi <- (dnorm(z/sigma)/(pnorm(z/sigma))) 
    dldv <- -(1/nu)+(z/nu^2)+(sigma/nu^2)*pphi
    dldv
                                    },
               d2ldv2 = function() {
       z <- y-mu-((sigma^2)/nu)
    pphi <- (dnorm(z/sigma)/(pnorm(z/sigma))) 
    dldv <- -(1/nu)+(z/nu^2)+(sigma/nu^2)*pphi
  d2ldv2 <- -dldv*dldv
                                    },
              d2ldmdd = function() rep(0,length(y)),
              d2ldmdv = function() rep(0,length(y)),
              d2ldddv = function() rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,nu,...) {
G.dev.incr <- -2*dexGAUS(y,mu=mu,sigma=sigma,nu=nu,log=TRUE)
                                                    }, 
             rqres = expression(rqres(pfun="pexGAUS", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
        mu.initial = expression( mu <- y+0.00001), 
     sigma.initial = expression( sigma <- rep(sd(y),length(y))), 
        nu.initial = expression( nu <- rep(sd(y), length(y))), 
          mu.valid = function(mu) TRUE , 
       sigma.valid = function(sigma)  all(sigma > 0),
          nu.valid = function(nu) all(nu > 0), 
           y.valid = function(y) TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
# the formula we use here is not robust for large sigma or small nu
# for small nu's it should go to normal I am not sure what should happent for 
# large sigmas
# the ratio sigma/nu is important
# i.e.
#> dexGAUS(1, mu=5, sigma=1, nu=.001)
#[1] 0
#> dexGAUS(1, mu=0, sigma=60, nu=.001)
#[1] 0
dexGAUS<-function(y, mu=5, sigma=1, nu=1, log=FALSE)
  { 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
 #  if (any(nu <= 1) )  stop(paste("nu must be greater than 1 ", "\n", "")) 
    ly <- length(y)       
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
     z <- y-mu-((sigma^2)/nu)
#cat("the zeta",z,"\n")
logfy <- ifelse(nu>0.001*sigma,  
          -log(nu)-(z+(sigma^2/(2*nu)))/nu+log(pnorm(z/sigma)),
           dnorm(y, mean=mu, sd=sigma, log=TRUE)
               )
#logfy <- ifelse(nu > 0.03, logfy, dnorm(y,mean=mu,sd=sigma, log=TRUE))
#logfy <- -log(nu)+((mu-y)/nu)+(sigma^2/(2*nu^2))+log(pnorm(((y-mu)/sigma)-sigma/nu))
  if(log==FALSE) fy <- exp(logfy) else fy <- logfy
  fy
  }
#----------------------------------------------------------------------------------------
# the some problem as above where the 
pexGAUS<-function(q, mu = 5, sigma = 1,  nu=1,  lower.tail = TRUE, log.p = FALSE)
  { 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(nu <= 0) )  stop(paste("nu must be greater than 0 ", "\n", "")) 
    ly <- length(q)       
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
 index <- seq(along=q)
     z <- q-mu-((sigma^2)/nu)
   cdf <- ifelse(nu>0.001*sigma,  
  pnorm((q-mu)/sigma)-pnorm(z/sigma)*exp(((mu+(sigma^2/nu))^2-(mu^2)-2*q*((sigma^2)/nu))/(2*sigma^2)),
          pnorm(q, mean=mu, sd=sigma)
               )
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#----------------------------------------------------------------------------------------
ppexGAUS<-function(q, mu = 5, sigma = 1,  nu=1,  lower.tail = TRUE, log.p = FALSE)
  { 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(nu <= 0) )  stop(paste("nu must be greater than 0 ", "\n", "")) 
    ly <- length(q)       
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
 index <- seq(along=q)
   cdf <-rep(0, ly)
 for (i in index)
          {
        cdf[i] <- integrate(function(x) 
                 dexGAUS(x, mu = mu[i], sigma = sigma[i], nu=nu[i], 
                      log=log.p), -Inf, q[i] )$value
          }    
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#----------------------------------------------------------------------------------------
# this is use golded section and is slower than the next version which use uniroot
#qexGAUS <- function(p, mu = 0, sigma = 1, nu = 1,  lower.tail = TRUE, 
#                    log.p = FALSE, lower.limit = mu-10*sqrt(sigma^2+nu^2),
#              upper.limit = mu+10*sqrt(sigma^2+nu^2))
#  { 
#    #---Golded section search--------------------------------------------
#       find.q.from.p <- function(p, mu, sigma, nu,
#                                 lower= lower.limit, 
#                                 upper = upper.limit)
#            {   
#               usemode <- function(q,p)
#                    { 
#                      np <- pexGAUS(q , mu = mu, sigma = sigma, nu = nu  )
#                      fun <- (np-p)^2
#                      fun
#                    }      
#            tol <-  0.000001
#              r <- 0.61803399  
#              b <- r*lower + (1-r)*upper 
#             lo <- lower
#             up <- upper 
#             w1 <- TRUE 
#           val1 <- if(up-b > b-lo) b   else b-(1-r)*(b-lo)
#           val2 <- if(up-b > b-lo) b+(1-r)*(up-b) else b 
#             f1 <- usemode(val1,p)
#             f2 <- usemode(val2,p)
#            while(w1) 
#                 { if(f2 < f1) { lo <- val1 
#                               val1 <- val2 
#                               val2 <- r*val1+(1-r)*up
#                                 f1 <- f2 
#                                 f2 <- usemode(val2,p) 
#                               } 
#                   else        { up <- val2 
#                               val2 <- val1
#                               val1 <- r*val2+(1-r)*lo
#                                 f2 <- f1 
#                                 f1 <- usemode(val1,p)
#                               }
#                   w1 <- abs(up-lo) >  tol*(abs(val1)+abs(val2))
#                 }             
#            q <- if(f1<f2) val1 else val2
#            q
#            }
#     #-----------------------------------------------------------------
#   if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
#    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
#    if (lower.limit >= upper.limit) stop(paste("the lower limit for the golden search is greated that the upper", "\n", ""))      
#    if (log.p==TRUE) p <- exp(p) else p <- p
#    if (lower.tail==TRUE) p <- p else p <- 1-p
#    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
#         lp <- length(p)                                                                    
#      sigma <- rep(sigma, length = lp)
#         mu <- rep(mu, length = lp)
#         nu <- rep(nu, length=lp)
#      upper <- rep(upper.limit, length = lp )
#      lower <- rep(lower.limit, length = lp )
#          q <- rep(0,lp)    
#         for (i in seq(along=p))
#          {
#          q[i] <- find.q.from.p(p[i], mu = mu[i], sigma = sigma[i], nu=nu[i],
#                                upper = upper[i], 
#                                lower = lower[i])
#          if (q[i]>=upper[i]) warning("q is at the upper limit, increase the upper.limit")
#          if (q[i]<=lower[i]) warning("q is at the lower limit, decrease the lower.limit")
#          }                                                                               
#    q
#   }
#----------------------------------------------------------------------------------------
qexGAUS <- function(p, mu = 5, sigma = 1, nu = 1,  lower.tail = TRUE, 
                    log.p = FALSE, lower.limit = mu-10*sqrt(sigma^2+nu^2),
              upper.limit = mu+10*sqrt(sigma^2+nu^2))
  { 
    usemode <- function(y)
                    { 
                      np <- pexGAUS(y , mu = mu[i], sigma = sigma[i], nu = nu[i]  )
                      fun <- (np-p[i])
                      fun
                    }     
     #-----------------------------------------------------------------
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(lower.limit >= upper.limit)) stop(paste("the lower limit is greated that the upper", "\n", ""))      
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
         lp <- length(p)                                                                    
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length=lp)
      upper <- rep(upper.limit, length = lp )
      lower <- rep(lower.limit, length = lp )
          q <- rep(0,lp)    
         for (i in seq(along=p))
          {
          interval <- c(lower[i],upper[i])
           q[i] <- if (usemode(interval[1]) * usemode(interval[2]) > 0) 
                          0
                    else uniroot(usemode, interval)$root
          if (q[i]>=upper[i]) warning("q is at the upper limit, increase the upper.limit")
          if (q[i]<=lower[i]) warning("q is at the lower limit, decrease the lower.limit")
          }                                                                               
    q
   }

#----------------------------------------------------------------------------------------
rexGAUS <- function(n, mu=5, sigma=1, nu=1, ...)
  { 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qexGAUS(p,mu=mu,sigma=sigma, nu=nu, ...)
    r
  }
#----------------------------------------------------------------------------------------
