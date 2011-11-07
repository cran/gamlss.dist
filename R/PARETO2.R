######################################
##### PARETO TYPE 2 DISTRIBUTION #####
######################################
#-------------------------------------------------------------------------------
# Robert Rigby, Fiona McElduf, Mikis Stasinopoulos, Vlasios Voudouris 
#################################################################################
#Gamlss Family Function
PARETO2 <- function (mu.link = "log", sigma.link = "log") 
{
    mstats <- checklink("mu.link", "Pareto Type 2", substitute(mu.link), 
        c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Pareto Type 2", substitute(sigma.link), 
        c("inverse", "log", "identity", "own"))
    structure(
       list(family = c("PARETO2", "Pareto Type 2"), 
        parameters = list(mu = TRUE, sigma = TRUE), 
             nopar = 2, 
              type = "Continuous", 
           mu.link = as.character(substitute(mu.link)), 
        sigma.link = as.character(substitute(sigma.link)), 
        mu.linkfun = mstats$linkfun, 
     sigma.linkfun = dstats$linkfun, 
        mu.linkinv = mstats$linkinv, 
     sigma.linkinv = dstats$linkinv, 
             mu.dr = mstats$mu.eta, 
          sigma.dr = dstats$mu.eta, 
              dldm = function(y, mu, sigma) #(1/(mu*sigma))-(1+(1/sigma))/(y+mu),
              {
               dldm1 <- (1/(mu*sigma))-(1+(1/sigma))*(1/(y+mu)) # bob's
               #dldm2 <- (y-mu*sigma)/(y*mu*sigma+mu^2*sigma) # fiona's 
                #   z <- (y/mu)
                #dldm3 <- -1/mu + (1/mu)*(1+(1/sigma))*z/(1+z) # gb2
               dldm1
                },      
            d2ldm2 = function(y, mu, sigma)
              {
              #d2ldm21 <- -(1/sigma)/((1+(1/sigma)+1)*mu^2) # gb2      
              #  dldm <-  (1/(mu*sigma))-(1+(1/sigma))/(y+mu) (1/(mu*sigma))-1/(y+mu)-(1/(sigma*(y+mu)))
              #d2ldm22 <- -(1/((mu^2)*sigma)) + (1+(1/sigma)/((y+mu)^2)) # bob's
              d2ldm23 <- -(1/(mu^2*(1+2*sigma))) # fiona's
            # browser() #1 and 3 similar
              d2ldm23
              },
              dldd = function(y, mu, sigma)  
              {
        #     z <- (y/mu)
        #   dldd1 <- -log(1+z)-digamma(1/sigma)+digamma(1+(1/sigma)) # bg2
        #   dldd     
           dldd2 <- -(1/sigma)+(log(y+mu)/(sigma^2))-(log(mu)/sigma^2) # bob's
        #   dldd3 <- (-sigma-log(mu)+log(y+mu))/(sigma^2) # fiona's
        #   cat(dldd[1], "\n")
         #   browser() #all different
            dldd2
            },
            d2ldd2 = function(y, mu, sigma)      
             {
              #dldd<- -(1/sigma)+(log(y+mu)/sigma^2)-(log(mu)/sigma^2)
             #d2ldd2<- -dldd^2
            #  cat(d2ldd2[1], "\n")
            #(1/sigma^2)+(2*log(mu))/(sigma^3)- (2*log(y+mu))/(sigma^3)
            # d2ldd21 <- -(1/sigma)+(log(y+mu)/(sigma^2))-(log(mu)/sigma^2) # bob's
            # d2ldd2 <- trigamma(1+(1/sigma)) - trigamma(1/sigma)
            d2ldd22 <- -(1/sigma^2) 
               #dldd2 <- -(1/sigma)+(log(y+mu)/(sigma^2))-(log(mu)/sigma^2)
               #d2ldd2 <- -dldd2^2
             d2ldd22  
              },
        d2ldmdd = function(y, mu, sigma) #
        {
          d2ldmdd <- -(1/(mu*sigma+mu*sigma^2))
       d2ldmdd 
        },
        G.dev.incr = function(y, mu, sigma, ...) -2 * 
            dPARETO2(y, mu, sigma, log = TRUE), 
        rqres = expression(rqres(pfun = "pPARETO2", 
            type = "Continuous", y = y, mu = mu, sigma = sigma)), 
        mu.initial = expression({mu <- (y + mean(y))/2}), 
        sigma.initial = expression({sigma <- rep(sd(y), length(y))}), 
        mu.valid = function(mu) all(mu > 0), 
        sigma.valid = function(sigma) all(sigma > 0), 
        y.valid = function(y) TRUE), 
        class = c("gamlss.family", "family"))
}
#-------------------------------------------------------------------------------
#Probability distribution function
dPARETO2 <- function(x, mu = 1, sigma = 0.5, log = FALSE)
{
   if (any(mu < 0)) stop(paste("mu must be positive", "\n", "")) 
   if (any(sigma <= 0))   stop(paste("sigma must be positive", "\n", ""))  
   if (any(x < 0)) stop(paste("x must be greater than 0", "\n", ""))
   lfy <- -log(sigma) + (1/sigma)*log(mu) - ((1/sigma)+1)*log(x+mu)
   if (log == FALSE) fy <- exp(lfy) else fy <- lfy
   fy
}
#--------------------------------------------------------------------------------
#Cumulative density function
pPARETO2 <- function(q, mu = 1, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
{
   if (any(mu <= 0)) stop(paste("mu must be positive", "\n", "")) 
   if (any(sigma <= 0)) stop(paste("tau must be positive", "\n", ""))           
   if (any(q < 0)) stop(paste("q must be be greater than 0", "\n", ""))   
   cdf <- 1 - ((mu/(mu+q))^(1/sigma))
   if (lower.tail == TRUE) cdf <- cdf  
   else cdf <- 1 - cdf
   if (log.p == FALSE) cdf <- cdf
   else cdf < - log(cdf)
   cdf
}   
#-------------------------------
#Quantile-inverse cdf  
qPARETO2 <- function(p, mu = 1, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
#    if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
   # w <- qf(p,2,2/sigma)   
   #  q1 <- mu*(((sigma)*w))
    q <-  mu*((1-p)^(-sigma)-1)    
    q
}
#--------------------------------------------------------------------------------
#Random generation 
rPARETO2 <- function(n, mu = 1, sigma = 0.5)
{
   if (any(mu <= 0)) stop(paste("mu must be positive", "\n", "")) 
   if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
   if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))  
   n <- ceiling(n)
   p <- runif(n)
   r <- qPARETO2(p, mu = mu, sigma = sigma)
   r 
}
#-------------------------------------------------------------------------------
