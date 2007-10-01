# Wednesday, January 21, 2004 at 18:25 MS BR
# the first derivatives squares have been used here
#  I need the variance of the function so I can put limits on the q function 
SHASH <- function (mu.link="identity", sigma.link="log", nu.link ="log", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Sinh-Arcsinh", substitute(mu.link), 
                           c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Sinh-Arcsinh", substitute(sigma.link), 
                           c("inverse", "log", "identity"))
    vstats <- checklink(   "nu.link", "Sinh-Arcsinh", substitute(nu.link),    
                           c("1/nu^2", "log", "identity"))
    tstats <- checklink(  "tau.link", "Sinh-Arcsinh", substitute(tau.link),   
                           c("1/tau^2", "log", "identity")) 
    structure(
          list(family = c("SHASH", "Sinh-Arcsinh"),
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
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm
                      },
   d2ldm2 = function(){
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm
 d2ldm2 <- -dldm*dldm
                      },     
   dldd = function() {  
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma                     
                      } ,
   d2ldd2 = function(){
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma    
 d2ldd2 <- -dldd*dldd
                      },   
     dldv = function() { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
      dldv <- dldr*drdv+dldc*dcdv
                        } ,
    d2ldv2 = function() { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
      dldv <- dldr*drdv+dldc*dcdv
      d2ldv2 <-  -dldv*dldv             
                        },
      dldt = function() {
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
      dldt <- dldr*drdt+dldc*dcdt
                        } ,
      d2ldt2 = function() { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
      dldt <- dldr*drdt+dldc*dcdt
      d2ldt2 <-   -dldt*dldt   
                            } ,
       d2ldmdd = function()## ok
               {
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm     
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma           
d2ldmdd <- -(dldm*dldd)
               },
       d2ldmdv = function()# OK
               { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm                   
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
   dldv <- dldr*drdv+dldc*dcdv
d2ldmdv <- -(dldm*dldv)
               },
       d2ldmdt = function() #ok
               {
          z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm     
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt                          
d2ldmdt <- -(dldm*dldt)
               },
       d2ldddv = function() #ok
               {               
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma       
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
   dldv <- dldr*drdv+dldc*dcdv
d2ldddv <- -(dldd*dldv)
               },
       d2ldddt = function() #ok
               {
     z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma   
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt   
d2ldddt <- -(dldd*dldt)  
               },
       d2ldvdt = function() #ok
               { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
   dldv <- dldr*drdv+dldc*dcdv                         
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt             
d2ldvdt <- -(dldv*dldt)  
               },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) -2*dSHASH(y,mu,sigma,nu,tau,log=TRUE),                 
         rqres = expression(   
                   rqres(pfun="pSHASH", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)
                           ),
    mu.initial = expression(mu <- (y+mean(y))/2),# rep(mean(y),length(y)) 
 sigma.initial = expression(sigma<- rep(sd(y)/5, length(y))),
    nu.initial = expression(nu <- rep(.5, length(y))), 
   tau.initial = expression(tau <-rep(.5, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dSHASH <- function(y, mu = 0, sigma = 1, nu = 1, tau = .5, log = FALSE)
 {
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
 loglik <- -log(sigma)-log(2*pi)/2-log(1+(z^2))/2+log(c)-(r^2)/2
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pSHASH <- function(q, mu = 0, sigma = 1, nu = 1, tau = .5, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
         if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          
      z <- (q-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
#     c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      p <- pNO(r)
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#-----------------------------------------------------------------  
qSHASH <-  function(p, mu=0, sigma=1, nu=1, tau=.5, lower.tail = TRUE, log.p = FALSE,
                     lower.limit = mu-10*(sigma/(nu*tau)), # this is completly wrong
                     upper.limit = mu+10*(sigma/(nu*tau)) )
 {   
     #---Golded section search--------------------------------------------
       find.q.from.p <- function(p, mu, sigma, nu, tau,
                                 lower= lower.limit, 
                                 upper = upper.limit)
            {   
               usemode <- function(q,p)
                    { 
                      np <- pSHASH(q , mu = mu, sigma = sigma, nu=nu, tau=tau  )
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
    #-----
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
         lp <-  pmax.int(length(p), length(mu), length(sigma), length(nu), length(tau))
          p <- rep(p, length = lp)
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length = lp)
        tau <- rep(tau, length = lp)
      upper <- rep(upper.limit, length = lp )
      lower <- rep(lower.limit, length = lp )
          q <- rep(0,lp)    
         for (i in seq(along=p))
          {
          q[i] <- find.q.from.p(p[i], mu = mu[i], sigma = sigma[i], nu=nu[i], tau=tau[i], 
                                upper = upper[i], 
                                lower = lower[i])
          if (q[i]>=upper[i]) warning("q is at the upper limit, increase the upper.limit")
          if (q[i]<=lower[i]) warning("q is at the lower limit, decrease the lower.limit")
          }                                                                               
    q
 }
#-----------------------------------------------------------------  
rSHASH <- function(n, mu=0, sigma=1, nu=1, tau=.5)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qSHASH(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
   r
  }
#-----------------------------------------------------------------  
