library(gamlss.dist)

## checking distributons
# two parameters 
# Continuous
#------------------------------------------------------------------------------------------
# Normal
#plot(function(y) dNO(y, mu=10 ,sigma=2), 0, 20)
#plot(function(y) pNO(y, mu=10 ,sigma=2), 0, 20)
#plot(function(y) qNO(y, mu=10 ,sigma=2), 0, 1)
#plot(function(y) qNO(y, mu=10 ,sigma=2,lower.tail=F), 0, 1)
#pNO(2,mu=10,sigma=2,lower.tail=T)+pNO(2,mu=10,sigma=2,lower.tail=F)
#pr <- pNO(2 ,mu=10,sigma=2)
#qq <- qNO(pr,mu=10,sigma=2)
#if(abs(qq-2)>0.0001) stop("error in NO") else cat("NO OK \n")
#aa<-integrate(function(x) dNO(x), -Inf, 1 )$value-pNO(1)
#if(abs(aa)>0.0001) stop("error in NO") else cat("NO OK \n")
#hist(rNO(100,mu=10,sigma=2))

#----------------------------------------------------------------------------------------
#-----WEI3  Weibull
integrate(function(x) dWEI3(x), 0, Inf)
plot(function(y) dWEI3(y, mu=10 ,sigma=2), 0, 30)
plot(function(y) pWEI3(y, mu=10 ,sigma=2), 0, 30)
plot(function(y) qWEI3(y, mu=10 ,sigma=2), 0, 1)
plot(function(y) qWEI3(y, mu=10 ,sigma=2, lower.tail=F), 0, 1)
pWEI3(7,mu=10,sigma=2,lower.tail=T)+pWEI3(7,mu=10,sigma=2,lower.tail=F)
pr <- pWEI3(10 ,mu=10,sigma=2)
qq <- qWEI3(pr ,mu=10,sigma=2)
if(abs(qq-10)>0.0001) stop("error in WEI3") else cat("WEI3 OK \n") 
aa<-integrate(function(x) dWEI3(x), 0, 2)$value-pWEI3(2)
if(abs(aa)>0.0001) stop("error in WEI3") else cat("WEI3 OK \n") 
hist(rWEI3(100,mu=10,sigma=2))
#----------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#--- Reverse Generalised extreme 
plot(function(y) dRGE(y, mu=0 ,sigma=1, nu=1), -1, 6)
plot(function(y) dRGE(y, mu=0 ,sigma=1, nu=.25), -3.999, 5)# 
plot(function(y) pRGE(y, mu=0 ,sigma=1, nu=1), -1, 6)

plot(function(y) qRGE(y, mu=0 ,sigma=1, nu=1), 0, 1)
plot(function(y) qRGE(y, mu=0 ,sigma=1, nu=1, lower.tail=F), 0, 1)


pRGE(1,mu=1,sigma=1, nu=1, lower.tail=T)+pRGE(1,mu=1,sigma=1, nu=1,lower.tail=F)
pr <- pRGE(3 ,mu=1,sigma=1, nu=1)
qq <- qRGE(pr ,mu=1,sigma=1, nu=1)
if(abs(qq-3)>0.0001) stop("error in RGE")  else cat("RGE OK \n") 

aa<-integrate(function(x) dRGE(x, mu=0, sigma=1, nu=1), -1, 2)$value-pRGE(2, mu=0, sigma=1, nu=1)
if(abs(aa)>0.0001) stop("error in RGE")  else cat("RGE OK \n") 
hist(rRGE(100,mu=0,sigma=1, nu=1))
#----------------------------------------------------------------------------------------

#----BEta
#plot(function(y) dBE(y, mu=.1 ,sigma=.5), 0.001, .999)
#plot(function(y) pBE(y, mu=.1 ,sigma=.5), 0.001, 0.999)
#plot(function(y) qBE(y, mu=.1 ,sigma=.5), 0.001, 0.999)
#plot(function(y) qBE(y, mu=.1 ,sigma=.5, lower.tail=F), 0.001, .999)
#pBE(.8, mu =.1, sigma=.9 ,lower.tail=T)+pBE(.8,mu=.1,sigma=.9,lower.tail=F)
#pr <- pBE(.8 ,mu=.1,sigma=.9)
#qq <- qBE(pr,mu=.1,sigma=.9)
#if(abs(qq-.8)>0.0001) stop("error in BE") else cat("BE OK \n") 
#aa<-integrate(function(x) dBE(x), 0, .8)$value-pBE(.8)
#if(abs(aa)>0.0001) stop("error in BE") else cat("BE OK \n") 
#hist(rBE(100,mu=.1,sigma=.5))
#----------------------------------------------------------------------------------------
## three paramemeter Continuous
## GG
plot(function(y) dGG(y, mu=5 ,sigma=.1, nu=.25), 0, 10) # alsmost symmetrical
plot(function(y) dGG(y, mu=5 ,sigma=.1, nu=-10), 0, 10) # right skew
plot(function(y) dGG(y, mu=5 ,sigma=.1, nu=10), 0, 10) # left skew
curve(dGG(y=x, mu=5 ,sigma=.5, nu=0), 0, 10) # log normal
plot(function(y) pGG(y, mu=5 ,sigma=.1, nu=10), 0, 10)
plot(function(y) qGG(y, mu=5 ,sigma=.1, nu=10), 0.001, .999)
plot(function(y) qGG(y, mu=5 ,sigma=.1, nu=10, lower.tail=FALSE), 0.001, .999)
plot(function(y) qGG(y, mu=5 ,sigma=.1, nu=0, lower.tail=TRUE), 0.001, .999)
plot(function(y) qGG(y, mu=5 ,sigma=.1, nu=0, lower.tail=FALSE), 0.001, .999)
pr <- pGG(4 ,mu=5,sigma=.5, nu=1)
qq <- qGG(pr,mu=5,sigma=.5, nu=1)
if(abs(qq-4)>0.0001) stop("error in   GG") else cat("GG OK \n") 
aa<-integrate(function(x) dGG(x), 0, 4 )$value-pGG(4)
if(abs(aa)>0.0001) stop("error in GG") else cat("GG OK \n") 
hist(rGG(100,mu=5,sigma=.1, nu=10))

#----------------------------------------------------------------------------------------
## GIG
plot(function(y) dGIG(y, mu=5 ,sigma=.1, nu=.25), 0, 10) #
#plot(function(y) dGIG(y, mu=5 ,sigma=.1, nu=-10), 0, 10) # 
#plot(function(y) dGIG(y, mu=5 ,sigma=.1, nu=10), 0, 10) # 
curve(dGIG(y=x, mu=5 ,sigma=1, nu=1), 0, 20)

plot(function(y) pGIG(y, mu=5 ,sigma=1, nu=1), 0, 10)
plot(function(y) qGIG(y, mu=5 ,sigma=.1, nu=1), 0.001, .999)
plot(function(y) qGIG(y, mu=5 ,sigma=.1, nu=1, lower.tail=FALSE), 0.001, .999)


pr <- pGIG(4 ,mu=5, sigma=.5, nu=1)
qq <- qGIG(pr,mu=5, sigma=.5, nu=1)
if(abs(qq-4)>0.0001) stop("error in   GIG") else cat("GIG OK \n") 
aa<-integrate(function(x) dGIG(x), 0, 4 )$value-pGIG(4)
if(abs(aa)>0.0001) stop("error in GIG") else cat("GIG OK \n") 
hist(rGIG(100,mu=5,sigma=.1, nu=1))


#----ST3
plot(function(y) dST5(y, mu=0 ,sigma=1, nu=0, tau=1), -10, 10)
curve(dST5(y=x, mu=0 ,sigma=1, nu=0, tau=1), -10, 10)
plot(function(y) pST5(y, mu=0 ,sigma=1, nu=0, tau=1), -10, 10)
plot(function(y) qST5(y, mu=0 ,sigma=1, nu=0, tau=1), 0.001, .999)
plot(function(y) qST5(y, mu=0 ,sigma=1, nu=0, tau=1, lower.tail=FALSE), 0.001, .999)
pr <- pST5(2 ,mu=1,sigma=1, nu=0, tau=1)
qq <- qST5(pr,mu=1,sigma=1, nu=0, tau=1)
if(abs(qq-2)>0.0001) stop("error in ST5") else cat("ST5 OK \n") 

pr <- pST5(2 ,mu=1,sigma=1, nu=1, tau=1)
qq <- qST5(pr,mu=1,sigma=1, nu=1, tau=1)
if(abs(qq-2)>0.0001) stop("error in ST5") else cat("ST5 OK \n") 

pr <- pST5(200 ,mu=1,sigma=1, nu=1, tau=1)
qq <- qST5(pr,mu=1,sigma=1, nu=1, tau=1)
if(abs(qq-200)>0.0001) stop("error in ST5") else cat("ST5 OK \n") 

aa<-integrate(function(x) dST5(x), -Inf, 1 )$value-pST5(1)
if(abs(aa)>0.0001) stop("error in ST5") else cat("ST5 OK \n") 
hist(rST5(100,mu=0,sigma=1, nu=.5, tau=1))
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# four parameters Continuous
# SHASH
plot(function(y) dSHASH(y, mu=0, sigma=1, nu=1, tau=2), -4, 4)
curve(dSHASH(x, mu=0, sigma=1, nu=1, tau=.3), -4, 4)
# this was to check how the variance depends on nu and tau
#ff <-numeric(10)
#for (i in 1:10)
# {
# ii<-integrate(function(y) (y^2)*dSHASH(y, mu=0, sigma=1, nu=1, tau=i), -Inf, Inf)
# cat(ff[i]<-ii[[1]],"\n")
# }
plot(function(y) pSHASH(y, mu=0, sigma=1, nu=1, tau=2), -4, 4)
plot(function(y) qSHASH(y, mu=0, sigma=1, nu=1, tau=2), 0.01, 0.99)
#plot(function(y) qSEP(y, mu=0 ,sigma=1, nu=-1, tau=2, lower.tail=FALSE), 0.01, 0.99)
pSHASH(1,mu=0,sigma=1, nu=1,tau=2,lower.tail=T)+ pSHASH(1,mu=0,sigma=1, nu=1,tau=2,lower.tail=F)
pr <- pSHASH(1 ,mu=0,sigma=1,nu=1, tau=2)
qq <- qSHASH(pr  ,mu=0,sigma=1,nu=1, tau=2)
if(abs(qq-1)>0.0001) stop("error in SHASH") else cat("SHASH OK \n")
aa<-integrate(function(x) dSHASH(x, mu=0, sigma=1, nu=1,tau=2), -Inf, 1 )$value-pSHASH(1 ,mu=0, sigma=1, nu=1,tau=2)
if(abs(aa)>0.0001) stop("error in SHASH")
hist(rSEP(100, mu=0,sigma=1, nu=1,tau=2))
#----------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------- 
# discrete distributions
# 2 parameters
#-- PO
#plot(function(y) dPO(y, mu=10 ), from=0, to=20, n=20+1, type="h") # pdf
#pdf.plot(family=PO, mu=10, min=0, max=20, step=1)                 # pdf  
#plot(function(y) pPO(y, mu=10 ), from=0, to=20, n=20+1, type="h") # cdf
#plot(seq(from=0,to=30),pPO(seq(from=0,to=30), mu=10), type="h")   # cdf
#ppPO <- pPO(seq(from=0, to=40), mu=10 ) # get cum probabilities
#if(!isTRUE(all.equal(qqPO <- qPO(ppPO, mu=10),seq(from=0, to=40)))) 
#           warning("possible problem with PIG")
#pPO(5 ,mu = 5, lower.tail = T) + pPO(5 , mu = 5, lower.tail = F) # should be one
#plot(ppPO, qPO(ppPO, mu=10, lower.tail=TRUE ), type="h")
#plot(ppPO, qPO(ppPO, mu=10, lower.tail=FALSE ), type="h")
#plot(function(y) qPO(y, mu=10 ), from=0.01, to=.99, n=20+1, type="h") # I have to think about this
#plot(function(y) qPO(y, mu=10, lower.tail=F ), from=0.01, to=.99, n=20+1, type="h") 
#pr <- pPO(5, mu = 5)
#qq <- qPO(pr, mu = 5)
#if(abs(qq-5)>0.0001) stop("error in PO") else cat("PO OK \n")
#aa<-sum(dPO(y=c(0:5), mu=5))-pPO(5 , mu = 5)
#if(abs(aa)>0.0001) stop("error in PO")
#tN <- table(Ni <- rPO(1000, mu=5))
#r <- barplot(tN, col='lightblue')
#--
#------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# ZIP2 distribution
# where mu is the mean
# pdf plot
plot(function(y) dZIP2(y, mu=10, sigma = 0.1 ), from=0, to=20, n=20+1, type="h") # pdf
pdf.plot(family=ZIP2, mu=10, sigma=0.1, min=0, max=20, step=1)                 # pdf  
# cdf plot
PPP <- par(mfrow=c(2,1))
plot(function(y) pZIP2(y, mu=10, sigma =0.1 ), from=0, to=20, n=20+1, type="h") # cdf
cdf<-pZIP2(0:20, mu=10, sigma=0.1) 
sfun1  <- stepfun(1:20, cdf, f = 0)
plot(sfun1, xlim=c(0,20), main="cdf(x)")
par(PPP)
#plot(seq(from=0,to=20),pZIP2(seq(from=0,to=20), mu=10, sigma=0.1), type="h")   # cdf
# get cdf values and chech if the inverse cdf produce similar results
ppZIP2 <- pZIP2(seq(from=0, to=30), mu=10, sigma=.1 ) # get cum probabilities
if(!isTRUE(all.equal(qqZIP2 <- qZIP2(ppZIP2,mu=10, sigma=.1), seq(from=0, to=30)))) 
           warning("possible problem with ZIP")
# should sum to one
pZIP2(5 ,mu = 5, sigma = 0.1, lower.tail = TRUE) + pZIP2(5 , mu = 5, sigma=0.1, lower.tail = FALSE) # #
# inverse ccdf plots 
plot(ppZIP2, qZIP2(ppZIP2, mu=10, sigma=.1, lower.tail=TRUE), type="h")
plot(ppZIP2, qZIP2(ppZIP2, mu=10, sigma=.1, lower.tail=FALSE ), type="h")
# more checks
pr <- pZIP2(5,  mu = 5, sigma=.1)
qq <- qZIP2(pr, mu = 5, sigma=.1)
if(abs(qq-5)>0.0001) stop("error in ZIP2") else cat("ZIP2 OK \n")
# 
aa<-sum(dZIP2(y=c(0:5), mu=5, sigma=.1))-pZIP2(5 , mu = 5, sigma=.1)
if(abs(aa)>0.0001) stop("error in ZIP2")
## plotting bar chart and ecdf
tN <- table(Ni <- rZIP2(1000, mu=5, sigma=.1))
PPP <- par(mfrow=c(2,1))
r <- barplot(tN, col='lightblue')
plot(ecdf(Ni))
par(PPP)
## example of step function
x=1:20
cdf<-pZIP2(c(0,x), mu=10, sigma=0.1) 
sfun1  <- stepfun(x, cdf, f = 0)
plot(sfun1)
#----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# discrete 3 parameters
# the Sichel
plot(function(y) dSICHEL(y, mu=10, sigma = 2 , nu=1), from=0, to=100, n=100+1, type="h") # pdf
pdf.plot(family=SICHEL,     mu=10, sigma = 2,  nu=1, min=0, max=100, step=1)

plot(function(y) pSI(y,     mu=10, sigma = 2 , nu=1), from=0, to=100, n=100+1, type="h") # cdf
plot(seq(from=0,to=100),pSICHEL(seq(from=0,to=100), mu=10, sigma=2, nu=1), type="h")   # cdf

sd10 <- sum(dSICHEL(0:10, mu=10, sigma=3, nu=4))
p10  <- pSICHEL(10, mu=10, sigma=3, nu=4)
 if (abs(sd10-p10)>0.000001) warning("possible problem with p SICHEL")


ppSICHEL <- pSICHEL(seq(from=0, to=50), mu=5, sigma=1 , nu=1) # get cum probabilities
if(isTRUE(!all.equal(qqSICHEL <- qSICHEL(ppSICHEL,mu=5,  sigma=1, nu=1), seq(from=0, to=50)))) 
       warning("possible problem with SICHEL")

plot(function(y) qSICHEL(y, mu=5, sigma=1, nu=1), from=0.01, to=.99, n=50+1, type="h") # I have to think about this

pr <- pSICHEL(5, mu = 5, sigma=1, n=1)
qq <- qSICHEL(pr, mu = 5, sigma=1, n=1)
if(abs(qq-5)>0.0001) stop("error in SICHEL") else cat("SICHEL OK \n")
aa<-sum(dSICHEL(y=c(0:5), mu=5, sigma=1, nu=1))-pSICHEL(5 , mu = 5, sigma=1, nu=1)
if(abs(aa)>0.0001) stop("error in SICHEL")
tN <- table(Ni <- rSICHEL(200, mu=5, sigma=1, nu=1))
r <- barplot(tN, col='lightblue')
# m1<-histDistN(Ni, "SICHEL")

#------------------------------------------------------------------------------------------
#curve(dDEL(y=x, mu=10, nu=.1), 0, 40, type="h")
#-----------------------------------------------------------------------------------------
# DEL distribution
# 
# pdf plot
plot(function(y) dDEL(y, mu=1, sigma = 0.5 , nu=.5), from=0, to=10, n=10+1, type="h") # pdf
pdf.plot(family=DEL, mu=5, sigma=0.5, nu=0.5, min=0, max=15, step=1)             # pdf  
# cdf plot
PPP <- par(mfrow=c(2,1))
plot(function(y) pDEL(y, mu=5, sigma =0.5, nu=0.5 ), from=0, to=20, n=20+1, type="h") # cdf
cdf<-pDEL(0:20, mu=5, sigma=0.5, nu=0.5) 
sfun1  <- stepfun(1:20, cdf, f = 0)
plot(sfun1, xlim=c(0,20), main="cdf(x)")
par(PPP)
#plot(seq(from=0,to=20),pZIP2(seq(from=0,to=20), mu=10, sigma=0.1), type="h")   # cdf
# get cdf values and chech if the inverse cdf produce similar results
ppDEL <- pDEL(seq(from=0, to=30), mu=10, sigma=.5, nu=0.5 ) # get cum probabilities
if(!isTRUE(all.equal(qqDEL <- qDEL(ppDEL,mu=10, sigma=.5, nu=.5), seq(from=0, to=30)))) 
           warning("possible problem with ZIP")
# should sum to one
pZIP2(5 ,mu = 5, sigma = 0.1, lower.tail = TRUE) + pZIP2(5 , mu = 5, sigma=0.1, lower.tail = FALSE) # #
# inverse ccdf plots 
plot(ppZIP2, qZIP2(ppZIP2, mu=10, sigma=.1, lower.tail=TRUE), type="h")
plot(ppZIP2, qZIP2(ppZIP2, mu=10, sigma=.1, lower.tail=FALSE ), type="h")
# more checks
pr <- pZIP2(5,  mu = 5, sigma=.1)
qq <- qZIP2(pr, mu = 5, sigma=.1)
if(abs(qq-5)>0.0001) stop("error in ZIP2") else cat("ZIP2 OK \n")
# 
aa<-sum(dZIP2(y=c(0:5), mu=5, sigma=.1))-pZIP2(5 , mu = 5, sigma=.1)
if(abs(aa)>0.0001) stop("error in ZIP2")
## plotting bar chart and ecdf
tN <- table(Ni <- rZIP2(1000, mu=5, sigma=.1))
PPP <- par(mfrow=c(2,1))
r <- barplot(tN, col='lightblue')
plot(ecdf(Ni))
par(PPP)
## example of step function
x=1:20
cdf<-pZIP2(c(0,x), mu=10, sigma=0.1) 
sfun1  <- stepfun(x, cdf, f = 0)
plot(sfun1)
#----------------------------------------------------------------------------------------
