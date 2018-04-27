#----------------------------------------------------------------#
######## NONPARAMETRIC - EXAM 2: dWADE ##########################
####### Mai-Anh Dang | Student ID: 21608631 
#
# The code includes:
#
# (1) dWADE: to estimate the determinants on House Price
# Data: Anglin.gency.1996

library(dplyr)
library(ggplot2)
library(np)
library(magrittr)
#------------------------------------------------------------------#


#### (1) Density-weighted average derivative estimator ###################

## Acquire data
dat <- read.csv('./data/anglin.gencay.1996.csv')
dat['sell'] = log(dat['sell'])  # transform
dat['lot'] = log(dat['lot'])
dat['bdms'] = log(dat['bdms'])
dat['fb'] = log(dat['fb'])
dat['sty'] = log(day['sty'])
View(dat)

x = dat %>% select(-sell)
p = dat$sell

## Defining the derivative of Gaussian kernels ##

ker <- function(v){  ## v is a scalar
  K <- 1/(sqrt(2*pi))*exp(-0.5*v^2) ## Gaussian
  return(K)
}

dker <- function(v){
  dK <- 1/(sqrt(2*pi))*exp(-0.5*v^2)*(-v) ## derivative of Gaussian 
  return(dK)
}

mvker <- function(V){ # V is a k-vector
  mv <- sapply(V, ker) # apply ker for each element in vector mv
  Kmv <- prod(mv) # multivariate kernel use the product of kernel
  return(Kmv)
}

mvdker <- function(X_i, X_j, h){ # x_i and x_j is k-vector
  k = length(X_i)
  u = (X_i - X_j)/h # u is a k-vector 
  mvdker <- rep(1, k)  
  for (i in (1:k)){
    mvdker[i] = dker(u[i])*mvker(u)/ker(u[i])
  }
  return(unlist(mvdker)) # results is as k-vector 
}

##### (1.1) Calculate bw ####


library(kedd)

hn = h.ccv(x$lot, deriv.order = 1)$h

##### (1.2) Estimate dWADE ####

n = length(p)
d = length(x[1,])
beta.dwade = rep(0, d)

for (i in(1:n)){
  for (j in (1:n)){
    if (j != i){
      beta.dwade = beta.dwade + p[i] * mvdker(x[i,], x[j,], hn)*(1/(hn^(d+1)))*(-2/(n*(n-1)))
    }
  }
}

dwade = beta.dwade / beta.dwade[1]


##### (1.3) OLS ####

lm <- lm(dat$sell ~ dat$lot + dat$bdms + dat$fb + dat$sty
         + dat$drv + dat$rec + dat$ffin + dat$ghw + dat$ca + dat$gar + dat$reg-1)
lm2 <- lm(dat$sell ~ dat$lot + dat$bdms + dat$fb + dat$sty
          + dat$drv + dat$rec + dat$ffin + dat$ghw + dat$ca + dat$gar + dat$reg)

summary(lm)
beta.ols<- lm$coefficients
ols = beta.ols / beta.ols[1]

ols_intercept = lm2$coefficients
ols_intercept = ols_intercept / ols_intercept[2]
ols_intercept 

res = cbind(ols_intercept, c(0,ols), c(0,dwade))
colnames(res) = c('ols_intercept', 'ols_no_int', 'beta_dwade')
res
#write.csv(res, './results/dwade_ols.csv')