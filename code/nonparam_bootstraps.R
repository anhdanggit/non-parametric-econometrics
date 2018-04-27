#----------------------------------------------------------------#
######## NONPARAMETRIC - EXAM 2: Boostraps ##########################
####### Mai-Anh Dang | Student ID: 21608631 
#
# The code includes:
#
# (1) estimate the densities of GDP by Least-Squared Cross-Validation
# (2) compute the bootstraps CI with Asym. Pivotal Test Statistics
# (3) use the boostrap CI for hypothesis testing
# 
# Data: GDP
#------------------------------------------------------------------#

library(dplyr)
library(ggplot2)
library(np)
library(magrittr)
library(readxl)
library(np)

# Read Data
dat <- read_excel('./data/GDP.xlsx')

X05 <- log(dat$`2005`)
X16 <- log(dat$`2016`)

## Define Kernel function

ker <- function(v,kern){
  if (kern=="Ep"){
    if (abs(v) <= 1){
      K<-0.75*(1-v^2)} ## [Epanechnikov]
    else {
      K <- 0
    }
  }
  if (kern=="Ga"){
    K<-1/(sqrt(2*pi))*exp(-0.5*v^2) ## [Gaussian]
  }
  return(K)
}

## sum of kernel
sumkern <- function(x,x0,h,kern){
  S <- 0
  for (i in seq(1,length(x),by=1)){
    arg<-(x[i]-x0)/h
    S<-S+ker(arg,kern)
  }
  return(S)
}

## sum of kernel squared
sumkernsqr <- function(x,x0,h,kern){
  S <- 0
  for (i in seq(1,length(x),by=1)){
    arg<-(x[i]-x0)/h
    S<-S+ker(arg,kern)^2
  }
  return(S)
}

## (2.1) Calculate the bw ####

## Least-Squared Cross-validation bw
cvls_05 <-npudensbw(X05,bwmethod="cv.ls",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_cvls_05 <- cvls_05$bw
cvls_16 <-npudensbw(X16,bwmethod="cv.ls",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_cvls_16 <- cvls_16$bw

## undersmoothed

h05 = h_cvls_05 - 0.0001
h16 = h_cvls_16 - 0.0001

## (2.2) Estimate f.hat.05 and f.hat.16 ####

x <- seq(min(X05),max(X05),by=0.02) # grid points
f.hat.05 <- rep(0, length(x))
f.hat.16 <- rep(0, length(x))
n = length(X05)

# estimate f.hat for each point in grid
for (i in (1:length(x))){
  f.hat.05[i]=(1/(n*h05))*sumkern(X05,x[i],h05,"Ep")
}

for (i in (1:length(x))){
  f.hat.16[i]=(1/(n*h16))*sumkern(X16,x[i],h16,"Ep")
}


## (2.3) Estimate f.hat.05 CI 95% ####

B=500

CI.lb.05=rep(0,length(x)) # CI lower bound
CI.ub.05=rep(0,length(x)) # CI upper bound

for (i in seq(1,length(x),1)){
  
  #Bootstrap "monte-carlo"
  Tstat.bs=rep(0,B)
  for (b in (1:B)){
    set.seed(b)
    X.bs=sample(X05, length(X05),replace=TRUE) # resampling, take the sample X, resample with n
    h.bs=h05
    f.hat.bs=(1/(n*h.bs))*sumkern(X.bs,x[i],h05,"Ep")
    Tstat.bs[b]=f.hat.bs-f.hat.05[i] # bootstrap statistic
    #print(b)
  } # End of Bootstrap loop
  crit=quantile(Tstat.bs,c(0.025,0.975)) # for two tailed 95% CI
  CI.lb.05[i]=f.hat.05[i]-crit[2]
  CI.ub.05[i]=f.hat.05[i]-crit[1]
  
} #end of outer loop


## Pivot Tsat bootstrap X05####

B = 500

CI.lb.05.pivot=rep(0,length(x)) # CI lower bound
CI.ub.05.pivot=rep(0,length(x)) # CI upper bound


for (i in (1:length(x))){
  
  Tstat.bs.05.pivot=rep(0,B) # store Tstat
  Sigma.Tstat.05 = rep(0, B) # store V^(1/2)
  
  for (b in (1:B)){ # boot-straps B times
    
    set.seed(b)
    X.bs.05=sample(X05, length(X05),replace=TRUE) # resampling X, to get X*
    f.hat.bs.05=(1/(n*h05))*sumkern(X.bs.05,x[i],h05,"Ep") # compute f.hat.bs for X*
    
    Sigma.Tstat.05[b] = sqrt((1/(n*(h05^2)))*sumkernsqr(X.bs.05, x[i], h05, "Ep") - (1/n)*(f.hat.bs.05)^2) 
    Tstat.bs.05.pivot[b]=f.hat.bs.05-f.hat.05[i] # bootstrap T-statistic
  } # End of Bootstrap loop
} #end of outer loop

rel_pivot <- data.frame(Tstat = Tstat.bs.05.pivot, Sigma = Sigma.Tstat.05)
rel_pivot['pivot_stat'] = rel_pivot['Tstat']/rel_pivot['Sigma']
rel_pivot = rel_pivot[apply(rel_pivot, 1, function(x) all(is.finite(x))), ] # only take the finite results
crit.pivot.05 = quantile(rel_pivot$pivot_stat, c(0.025, 0.975)) # critical value by pivot

sd.f.x = rep(0, length(x))
for (i in (1:length(x))){
  sd.f.x[i] = sqrt((1/(n*(h05^2)))*sumkern(X05, x[i], h05, "Ep") - (1/n)*(f.hat.05)^2)
  CI.lb.05.pivot[i] = f.hat.05[i]-crit.pivot.05[2]*sd.f.x[i] # appro. student dist, symmetry
  CI.ub.05.pivot[i] = f.hat.05[i]+crit.pivot.05[2]*sd.f.x[i]
}

## (2.4) Estimate f.hat.16 CI 95% ####

B=500

CI.lb.16=rep(0,length(x)) # CI lower bound
CI.ub.16=rep(0,length(x)) # CI upper bound

for (i in seq(1,length(x),1)){
  
  #Bootstrap "monte-carlo"
  Tstat.bs.16=rep(0,B)
  for (b in (1:B)){
    set.seed(b)
    X.bs.16=sample(X16, length(X16),replace=TRUE) # resampling, take the sample X, resample with n
    h.bs=h16
    f.hat.bs.16=(1/(n*h.bs))*sumkern(X.bs.16,x[i],h16,"Ep")
    Tstat.bs.16[b]=f.hat.bs.16-f.hat.16[i] # bootstrap statistic
    #print(b)
  } # End of Bootstrap loop
  crit.16=quantile(Tstat.bs.16,c(0.025,.975)) # for two tailed 95% CI
  CI.lb.16[i]=f.hat.16[i]-crit.16[2]
  CI.ub.16[i]=f.hat.16[i]-crit.16[1]
  
} #end of outer loop

## Pivot Tsat bootstrap X16####

B = 500

CI.lb.16.pivot=rep(0,length(x)) # CI lower bound
CI.ub.16.pivot=rep(0,length(x)) # CI upper bound


for (i in seq(1,length(x),1)){
  
  Tstat.bs.16.pivot=rep(0,B) # store Tstat
  Sigma.Tstat.16 = rep(0, B) # store V^(1/2)
  
  for (b in (1:B)){ # boot-straps B times
    
    set.seed(b)
    X.bs.16=sample(X16, length(X16),replace=TRUE) # resampling X, to get X*
    f.hat.bs.16=(1/(n*h16))*sumkern(X.bs.16,x[i],h16,"Ep") # compute f.hat.bs for X*
    
    Sigma.Tstat.16[b] = sqrt((1/(n*(h16^2)))*sumkernsqr(X.bs.16, x[i], h16, "Ep") - (1/n)*(f.hat.bs.16)^2) 
    Tstat.bs.16.pivot[b]=f.hat.bs.16-f.hat.16[i] # bootstrap T-statistic
  }
} #end of outer loop

rel_pivot <- data.frame(Tstat = Tstat.bs.16.pivot, Sigma = Sigma.Tstat.16)
rel_pivot['pivot_stat'] = rel_pivot['Tstat']/rel_pivot['Sigma']
rel_pivot = rel_pivot[apply(rel_pivot, 1, function(x) all(is.finite(x))), ]
crit.pivot.16 = quantile(rel_pivot$pivot_stat, c(0.025, 0.975)) # critical value by pivot

sd.f.x = rep(0, length(x))
for (i in (1:length(x))){
  sd.f.x[i] = sqrt((1/(n*(h16^2)))*sumkern(X16, x[i], h16, "Ep") - (1/n)*(f.hat.16)^2)
  CI.lb.16.pivot[i] = f.hat.16[i]-crit.pivot.16[2]*sd.f.x[i] # appro. student dist, symmetry
  CI.ub.16.pivot[i] = f.hat.16[i]+crit.pivot.16[2]*sd.f.x[i]
}

#### Plot CI 95% ##########
dev.new()
rafalib::mypar(2,)

matplot(x,cbind(f.hat.05,CI.lb.05.pivot,CI.ub.05.pivot, CI.ub.05, CI.lb.05),
        type="l",col=c("black","red","red","blue", "blue"),lty=c(1,1,1,1,1),lwd = c(2,1,1,1,1),
        main = expression(paste("(a) 95% Confidence Intervals for the Density of GDP 2005")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log(GDP)")))
legend("topright", legend=c("t-stat", "Pivotal t-stat"), col=c("blue", "red"), lty=c(1,1))

matplot(x,cbind(f.hat.16, CI.lb.16.pivot,CI.ub.16.pivot, CI.ub.16, CI.lb.16),
        type="l",col=c("black","red","red","blue", "blue"),lty=c(1,1,1,1,1),lwd = c(2,1,1,1,1),
        main = expression(paste("(b) 95% Confidence Intervals for the Density of GDP 2016")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log(GDP)")))
legend("topright", legend=c("t-stat", "Pivotal t-stat"), col=c("blue", "red"), lty=c(1,1))


## (2.5) Testing distribution of X05 and X16 ####

#### X05 ~ Gamma (shape = 2, scale = 1) ####
dev.new()
rafalib::mypar(2,2)
testpdf01 = dgamma(x,scale = 1, shape = 2)
matplot(x,cbind(testpdf01,CI.lb.05.pivot, CI.ub.05.pivot),type="l",
        col=c("red","blue","blue"),lty=c(2,1,1), lwd = c(2,1,1),
        main = expression(paste("(i) H0: X05 ~ Gamma, scale = 1, shape = 2")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log(GDP)")))
legend("topright", legend=c("Null", "X05 CI"), col=c("red", "blue"), lty=c(2,1))

#### X05 ~ N (2,2) ####
testpdf02= dnorm(x, 2, sd = sqrt(2))
matplot(x,cbind(testpdf02,CI.lb.05.pivot, CI.ub.05.pivot),type="l",
        col=c("red","blue","blue"),lty=c(2,1,1), lwd = c(2,1,1),
        main = expression(paste("(ii) H0: X05 ~ N (2, 2)")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log(GDP)")))
legend("topright", legend=c("Null", "X05 CI"), col=c("red", "blue"), lty=c(2,1))

#### X16 ~ Gamma (shape = 2.5, scale = 0.5) ####
testpdf03 = dgamma(x,scale = 0.5, shape = 2.5)
matplot(x,cbind(testpdf03,CI.lb.16.pivot, CI.ub.16.pivot),type="l",
        col=c("red","blue","blue"),lty=c(2,1,1), lwd = c(2,1,1),
        main = expression(paste("(iii) H0: X16 ~ Gamma, scale = 0.5, shape = 2.5")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log(GDP)")))
legend("topright", legend=c("Null", "X16 CI"), col=c("red", "blue"), lty=c(2,1))

#### X16 ~ N (2,2) ####
testpdf04= dnorm(x, 3.75, sd = sqrt(2.75))
matplot(x,cbind(testpdf04,CI.lb.16.pivot, CI.ub.16.pivot),type="l",
        col=c("red","blue","blue"),lty=c(2,1,1), lwd = c(2,1,1),
        main = expression(paste("(iv) H0: X16 ~ N (3.75, 2.75)")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log(GDP)")))
legend("topright", legend=c("Null", "X16 CI"), col=c("red", "blue"), lty=c(2,1))

#### X05 = X016 ####
dev.new()
matplot(x,cbind(f.hat.16,CI.lb.05.pivot, CI.ub.05.pivot),type="l",
        col=c("red","blue","blue"),lty=c(2,1,1), lwd = c(2,1,1),
        main = expression(paste("(v) H0: X05 = X16")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log(GDP)")))
legend("topright", legend=c("X16", "X05 CI"), col=c("red", "blue"), lty=c(2,1))