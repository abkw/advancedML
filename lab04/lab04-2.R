library(kernlab)
library(AtmRay)
library("mvtnorm")


#The squared exponential function###############################################
SquaredExpKernel <- function(x1,x2,sigmaF=1,l=0.3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}



# Simulates nSim realizations (function) from a GP with mean m(x) and covariance K(x,x')
# over a grid of inputs (x)
SimGP <- function(m = muPrior,K,x,nSim,...){
  n <- length(x)
  if (is.numeric(m)) meanVector <- rep(0,n) else meanVector <- m(x)
  covMat <- K(x,x,...)
  f <- rmvnorm(nSim, mean = meanVector, sigma = covMat)
  return(f)
}


#####################Setting the parameters#####################################
sigmaF = 1
l = 0.3
xy = c(0.4, 0.719)
sigmaNoise = 0.1
muPrior = c(0,0)
xGrid <- seq(-1,1,length=20)
nSim = 1


fSim <- SimGP(m=muPrior, K=SquaredExpKernel, x=xGrid, nSim, sigmaF, l)
plot(xGrid, fSim[1,], type="p", ylim = c(-1,1))
if(nSim>1){
  for (i in 2:nSim) {
    lines(xGrid, fSim[i,], type="l")
  }
}


#Question1:  The posterior function ########################################################
posteriorGP = function(X, y, XStar, sigmaNoise, k){
  #The algorithm
  K = k(X, X)
  n = length(X)
  
  L = chol(K + sigmaNoise^2*diag(n))
  alpha = solve(t(L),solve(L,y))
  
  if (length(XStar)==1){
    KStar = XStar
  }
  else{
    KStar = k(X,XStar)  
  }
  
  FStar = t(KStar)*alpha
  v = solve(L,KStar)
  VFStar = k(XStar, XStar)-t(v)%*%v
  logP = -0.5*t(y)%*%alpha - sum(log(diag(L))) - (n/2)*log(2*pi)
  
  return(list('mean'=FStar, 'variance'= VFStar, 'logMargLik' = logP))
}


#Question2: evaluating on (0.4,0.719)###########################################
x = 0.4
y = 0.719
posteriorValue = posteriorGP(X = x,
                             y = y,
                             XStar = xGrid,
                             sigmaNoise = sigmaNoise,
                             k = SquaredExpKernel)
postMean = posteriorValue$mean
postVar = posteriorValue$variance

#Plotting the posterior mean
plotMean = function(mean,var){
  plot(xGrid, mean, type="p", ylab="f(x)", xlab="x", ylim = c(-2,2))
  
  #Plotting confidence band
  lines(xGrid, mean - 1.96*sqrt(diag(var)), col = "blue", lwd = 2)
  lines(xGrid, mean + 1.96*sqrt(diag(var)), col = "blue", lwd = 2)
}

plotMean(postMean,postVar)

#Question3: updating the posterior with observation (-0.6, -0.044)

x = -0.6
y = -0.044

posteriorValue2 = posteriorGP(X = x,
                              y = y,
                              XStar = xGrid,
                              sigmaNoise = sigmaNoise,
                              k = SquaredExpKernel)
postMean2 = posteriorValue2$mean
postVar2 = posteriorValue2$variance

par(mfrow = c(1,2))
plot(xGrid,postMean2)
plotMean(postMean2, postVar2)

#Question4: The posterior with five observations
x = c(-1.0,-0.6,-0.2,0.4,0.8)
y = c(0.768,-0.044,-0.940,0.719,-0.664)

posteriorValue3 = posteriorGP(X = x,
                              y = y,
                              XStar = xGrid,
                              sigmaNoise = sigmaNoise,
                              k = SquaredExpKernel)
postMean3 = posteriorValue3$mean
postVar3 = posteriorValue3$variance

plot(xGrid, rowSums(postMean3), type="l", ylab="f(x)", xlab="x", ylim = c(-2.5,2.5))

#Plotting confidence band
lines(xGrid, rowSums(postMean3) - 1.96*sqrt(diag(postVar3)), col = "blue", lwd = 2)
lines(xGrid, rowSums(postMean3) + 1.96*sqrt(diag(postVar3)), col = "blue", lwd = 2)






