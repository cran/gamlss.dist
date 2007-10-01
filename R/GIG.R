GIG <- function (mu.link="log", sigma.link="log", nu.link ="identity") 
{
    mstats <- checklink("mu.link", "GIG", substitute(mu.link), 
                         c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "GIG", substitute(sigma.link), #
                         c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "GIG",substitute(nu.link), 
                         c("1/nu^2", "log", "identity"))  
    
    structure(
          list(family = c("GIG", "Generalised Inverse Gaussian"),
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
       c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))  
    dldm <- -(nu/mu)+((c*y)/(2*sigma*mu^2))-1/(2*sigma*c*y)
    dldm
                                    },
        
           d2ldm2 = function() {
         c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))  
    d2ldm2 <- (nu-(c/sigma))/(mu^2)
    d2ldm2
                                    },

                 dldd = function() {
       c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))  
    dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)
    dldd <- 1/sigma*(nu-c/sigma+(1/(2*sigma))*((c*y/mu)+(mu/(c*y)))+dcdd*(sigma*nu/c-(1/2)*(y/mu-mu/(c^2*y))))    
    dldd
                                    },
               d2ldd2 = function() {
     #      this needs checking 
     #    c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))  
     # dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)    
    #d2cdd2<-(c*(2*nu+1)+dcdd*(2*sigma*nu-2*c-sigma))/(sigma*sigma)
    #d2ldcdd <- 1/sigma*((-1/sigma)+(1/(2*sigma))*((y/mu)-(mu/(c^2*y)))-dcdd*(((sigma*nu)/c^2)+mu/(c^3*mu)))
    #d2ldd2 <- 1/sigma*((nu+1)/sigma+dcdd*nu-(sigma/c)*d2cdd2+d2ldcdd*dcdd)
    #d2ldd2
    # at the momment have the first derivative square
        c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))  
     dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)
     dldd <- 1/sigma*(nu-c/sigma+(1/(2*sigma))*((c*y/mu)+(mu/(c*y)))+dcdd*(sigma*nu/c-(1/2)*(y/mu-mu/(c^2*y))))    
    d2ldd2 <- -dldd*dldd
    d2ldd2  
                                    },
                 dldv = function() {
       nd <- numeric.deriv(dGIG(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
     dldv <- as.vector(attr(nd, "gradient"))                
     dldv
                                    },
               d2ldv2 = function() {
        nd <- numeric.deriv(dGIG(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
      dldv <- as.vector(attr(nd, "gradient"))                
    d2ldv2 <- -dldv*dldv
    d2ldv2
                                    },
              d2ldmdd = function() {
        c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))  
     dldm <- -nu/mu+c*y/(2*sigma*mu^2)-1/(2*sigma*c*y)
     dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)
   # dcdd <- 1/c*(nu-1/2*sigma*(c*y/mu-mu/(c*y)))
     dldd <- 1/sigma*(nu-c/sigma+1/(2*sigma)*(c*y/mu+mu/(c*y))+dcdd*(sigma*nu/c-1/2*(y/mu-mu/(c^2*y))))
  d2ldmdd <- -dldm*dldd
  d2ldmdd
                        },
              d2ldmdv = function() { 
       c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))   
    dldm <- -nu/mu+c*y/(2*sigma*mu^2)-1/(2*sigma*c*y)
      nd <- numeric.deriv(dGIG(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
    dldv <- as.vector(attr(nd, "gradient"))  
 d2ldmdv <- -dldm*dldv      
 d2ldmdv        
                        },
              d2ldddv = function() {
       c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))   
    dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)
    dldd <- 1/sigma*(nu-c/sigma+1/(2*sigma)*(c*y/mu+mu/(c*y))+dcdd*(sigma*nu/c-1/2*(y/mu-mu/(c^2*y))))
    dldv <- as.vector(attr(nd, "gradient"))
 d2ldddv <- -dldv*dldd  
 d2ldddv
                        },
          G.dev.incr  = function(y,mu,sigma,nu,...) 
                     -2*dGIG(y,mu=mu,sigma=sigma,nu=nu,log=TRUE), 
             rqres = expression(rqres(pfun="pGIG", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
        mu.initial = expression( mu <- (y+mean(y))/2), 
     sigma.initial = expression( sigma <- rep(1, length(y))), # sd(y)/(mean(y))^1.5 ), 
        nu.initial = expression( nu <- rep(1, length(y))), 
          mu.valid = function(mu) TRUE , 
       sigma.valid = function(sigma)  all(sigma > 0),
          nu.valid = function(nu) TRUE , 
           y.valid = function(y) all(y>0)
          ),
            class = c("gamlss.family","family"))
}
#--------------------------------------------------------------
dGIG <- function(y, mu=1, sigma=1, nu=1,  log = FALSE)
 {
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(y < 0))  stop(paste("y must be positive", "\n", "")) 
               c <- exp(log(besselK(1/sigma,nu+1))-log(besselK(1/sigma,nu)))  
          loglik <- nu*log(c)-nu*log(mu)+(nu-1)*log(y)-log(2)-log(besselK(1/sigma,nu))-1/(2*sigma)*(c*y/mu+mu/(c*y))
          if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
          ft
  }    
#--------------------------------------------------------------  
pGIG <- function(q, mu=1, sigma=1, nu=1,  lower.tail = TRUE, log.p = FALSE)
 {  
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("q must be positive", "\n", ""))  
         lq <- length(q)       
      sigma <- rep(sigma, length = lq)
         mu <- rep(mu, length = lq)   
         nu <- rep(nu, length = lq) 
        cdf <-rep(0, lq)
       for (i in 1:lq)
          {
        cdf[i] <- integrate(function(x) 
                 dGIG(x, mu = mu[i], sigma = sigma[i], nu = nu[i], log=log.p), 0.001, q[i] )$value
          }    
        cdf
 }
#--------------------------------------------------------------
qGIG <- function(p,  mu=1, sigma=1, nu=1, lower.tail = TRUE, log.p = FALSE, 
                lower.limit = 0,
                upper.limit = mu+10*sqrt(sigma^2*mu^3))
 { 
 #------
 find.q.from.p <- function(p, mu, sigma, 
                                 lower= lower.limit, 
                                 upper = upper.limit)
            {   
               usemode <- function(q,p)
                    { 
                       np <- pGIG(q , mu = mu, sigma = sigma  )
                      fun <- (np-p)^2
                      fun
                    }      
            tol <-  0.000001
              r <- 0.61803399  
              b <- r*lower + (1-r)*upper 
             lo <- lower
             up <- upper 
             w1 <- TRUE 
           val1 <- if(up-b > b-lo) b   else b-(1-r)*(b-lo)
           val2 <- if(up-b > b-lo) b+(1-r)*(up-b) else b 
             f1 <- usemode(val1,p)
             f2 <- usemode(val2,p)
            while(w1) 
                 { if(f2 < f1) { lo <- val1 
                               val1 <- val2 
                               val2 <- r*val1+(1-r)*up
                                 f1 <- f2 
                                 f2 <- usemode(val2,p) 
                               } 
                   else        { up <- val2 
                               val2 <- val1
                               val1 <- r*val2+(1-r)*lo
                                 f2 <- f1 
                                 f1 <- usemode(val1,p)
                               }
                   w1 <- abs(up-lo) >  tol*(abs(val1)+abs(val2))
                 }             
            q <- if(f1<f2) val1 else val2
            q
            }
    #---------- 
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
    if (lower.limit >= upper.limit) stop(paste("the lower limit for the golden search is greated that the upper", "\n", ""))      
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
         lp <- length(p)                                                                    
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
      upper <- rep(upper.limit, length = lp )
      lower <- rep(lower.limit, length = lp )
          q <- rep(0,lp)    
         for (i in seq(along=p))
          {
          q[i] <- find.q.from.p(p[i], mu = mu[i], sigma = sigma[i], 
                                upper = upper[i], 
                                lower = lower[i])
          if (q[i]>=upper[i]) warning("q is at the upper limit, increase the upper.limit")
          if (q[i]<=lower[i]) warning("q is at the lower limit, decrease the lower.limit")
          }                                                                               
    q
 }
#--------------------------------------------------------------
rGIG <- function(n, mu=1, sigma=1, nu=1, ...)
  {
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qGIG(p,mu=mu,sigma=sigma,nu=nu, ...)
    r
  }
#--------------------------------------------------------------
