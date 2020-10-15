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


#####################Setting the parameters#####################################
sigmaF = 1
l = 0.3
sigmaNoise = 0.1
xGrid <- seq(-1,1,length=20)
nSim = 1


#Question1:  The posterior function ########################################################
posteriorGP = function(X, y, XStar, sigmaNoise, k){
  #The algorithm
  K = k(X, X)
  n = length(X)
  
  L = chol(K + sigmaNoise^2*diag(n))
  alpha = solve(t(L),solve(L,y))
  KStar = k(X,XStar)
  FStar = t(KStar)%*%alpha
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
plotMean = function(mean,var, title){
  plot(xGrid, mean, type="p", ylab="f(x)", xlab="x", ylim = c(-2,2), main = title)
  
  #Plotting confidence band
  lines(xGrid, mean - 1.96*sqrt(diag(var)), col = "blue", lwd = 2)
  lines(xGrid, mean + 1.96*sqrt(diag(var)), col = "blue", lwd = 2)
}
par(mfrow = c(1,1))
plotMean(postMean,postVar, title = 'x = 0.4, y = 0.719')

#Question3: updating the posterior with observation (-0.6, -0.044)

x = c(0.4,-0.6)
y = c(0.719,-0.044)

posteriorValue2 = posteriorGP(X = x,
                              y = y,
                              XStar = xGrid,
                              sigmaNoise = sigmaNoise,
                              k = SquaredExpKernel)
postMean2 = posteriorValue2$mean
postVar2 = posteriorValue2$variance

plotMean(postMean2, postVar2, title = 'x = (0.4,-0.6), y =(0.719,-0.044)')


#Question4: The posterior with five observations################################
x = c(0.4, -0.6, -1.0, -0.6, -0.2, 0.4, 0.8)
y = c(0.719, -0.044, 0.768, -0.044, -0.940, 0.719, -0.664)

posteriorValue3 = posteriorGP(X = x,
                              y = y,
                              XStar = xGrid,
                              sigmaNoise = sigmaNoise,
                              k = SquaredExpKernel)
postMean3 = posteriorValue3$mean
postVar3 = posteriorValue3$variance

plotMean(postMean3,postVar3, title = "All X and Y plot")



#Question5: Repeating with different hyperparameters.

sigmaF = 1
L = 1
posteriorValue4 = posteriorGP(X = x,
                              y = y,
                              XStar = xGrid,
                              sigmaNoise = sigmaNoise,
                              k = SquaredExpKernel)
postMean4 = posteriorValue4$mean
postVar4 = posteriorValue4$variance

for (i in 1:length(x)){
  plotMean(postMean4[,i],postVar4, title = paste('x = ',x[i], ',  y = ',y[i]))
}




