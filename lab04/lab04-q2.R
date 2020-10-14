data = read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv", header=TRUE, sep=";")

date = data$date
time = c(1:dim(data)[1])
temp = data$temp
day = rep(c(1:365), length(time)/365)


#Takinng every fifth observation

time = seq(from = 1, to = length(time), by = 5)
day = day[time]

scaledTime = (time-mean(time))/sd(time)
#Question1: familirize yourserl#################################################
library(kernlab)


#Defining the exponnetial kernel function

Matern32 <- function(sigmaf = 1, ell = 1) 
{
  rval <- function(x, y = NULL) {
    r = sqrt(crossprod(x-y));
    return(sigmaf^2*(1+sqrt(3)*r/ell)*exp(-sqrt(3)*r/ell))
  }
  class(rval) <- "kernel"
  return(rval)
} 

#Evaluating the function
MaternFunc = Matern32(sigmaf = 20, ell = 0.2) # MaternFunc is a kernel FUNCTION
MaternFunc(1,2)

#Computing the covarance matrix K
X = c(1,3,4)
Xstar = c(2,3,4)

K <- kernelMatrix(kernel = MaternFunc, x = X, y = Xstar)
K


#Question2: ####################################################################

#Letting sigma noise equal to residual variance of a quadratic regression fit using lm
modelData = data.frame(data$temp[time], scaledTime)
colnames(modelData) = c('temp', 'time')
polyFit <- lm(temp ~ time + I(time^2) + I(time^3), data = modelData)
sigmaNoise = sd(polyFit$residuals)
sigmaNoise


plot(scaledTime,data$temp[time], ylim = c(-25,30))

# Fit the GP with built in Square expontial kernel (called rbfdot in kernlab)

MaternFunc = Matern32(sigmaf = 20, ell = 0.2)
GPfit <- gausspr(scaledTime, data$temp[time], 
                 kernel = MaternFunc,
                 kpar = list(sigmaf = 20, ell = 0.2),
                 var = sigmaNoise^2)
meanPred <- predict(GPfit, scaledTime)
lines(scaledTime, meanPred, col="blue", lwd = 2)


#Question3: kernlab ############################################################
ell <- 0.5
SEkernel <- rbfdot(sigma = 1/(2*ell^2))

x<-scaledTime
xs<-scaledTime # XStar
n <- length(x)
Kss <- kernelMatrix(kernel = SEkernel, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = SEkernel, x = x, y = x)
Kxs <- kernelMatrix(kernel = SEkernel, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs) # Covariance matrix of fStar

# Probability intervals for fStar
lines(xs, meanPred - 1.96*sqrt(diag(Covf)), col = "green", lwd = 2)
lines(xs, meanPred + 1.96*sqrt(diag(Covf)), col = "green", lwd = 2)

# Prediction intervals for yStar
lines(xs, meanPred - 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "brown")
lines(xs, meanPred + 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "brown")

legend('bottomright',
       legend = c('Post Mean of First Model', 'fStar prob interval', 'fStar pred interval', 'Post Mean of Second Model'),
       col = c('blue', 'green', 'brown', 'red'),
       lty=1:2, 
       cex=0.8)



#Question4: ####################################################################

#Letting sigma noise equal to residual variance of a quadratic regression fit using lm
modelData = data.frame(data$temp[time], day)
colnames(modelData) = c('temp', 'day')
polyFit <- lm(temp ~ day + I(day^2) + I(day^3), data = modelData)
sigmaNoise = sd(polyFit$residuals)
sigmaNoise

# Fit the GP with built in Square expontial kernel (called rbfdot in kernlab)

MaternFunc = Matern32(sigmaf = 20, ell = 0.2)
GPfit <- gausspr(day,
                 data$temp[time], 
                 kernel = MaternFunc, 
                 kpar = list(sigmaf = 20, ell = 0.2),
                 var = sigmaNoise^2)
meanPred <- predict(GPfit, day)
lines(scaledTime, meanPred, col="red", lwd = 2)

#The first model has a smoother posterior mean than the second model, the second model seems not smooth.

#Question5: ####################################################################

GPKernel = function(l1,l2, sigmaF, d){
  
  rval <- function(x1, x2 = NULL) {
    return(sigmaF^2 * exp(-(2*(sin(pi*abs(x1-x2)/d))^2)/l1^2)*(exp(-0.5*abs(x1-x2)^2)/l2^2))
  }
  class(rval) <- "kernel"
  return(rval)
  
}

l1 = 1
l2 = 10
sigmaF = 20
d = 365/sd(time)

GPKernelFun = GPKernel(l1,l2,sigmaF, d)
GPKernelFun(1,2)

#Fitting the generlaized model
GPfit <- gausspr(scaledTime, data$temp[time], 
                 kernel = GPKernelFun,
                 kpar = list(sigmaF = sigmaF, l1 = l1, l2 = l2, d=d),
                 var = sigmaNoise^2)
meanPred <- predict(GPfit, scaledTime)
plot(time,data$temp[time], ylim = c(-25,30), main = 'Generalized Model')
lines(time, meanPred, col="blue", lwd = 2)
