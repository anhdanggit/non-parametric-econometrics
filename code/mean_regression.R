#----------------------------------------------------------------#
######## NONPARAMETRIC - EXAM 1: Mean Regression Functions #######
####### Mai-Anh Dang | Student ID: 21608631 
#
# The code includes:
#
# estimate log(income) = g(log(food expenditure))
#
# (1) Local Constant Estimator
# (2) Local Linear Estimator
# (3) OLS
# 
# Data: Engel.dta
#------------------------------------------------------------------#



##### Obtain Data ####

library(readstata13)
engel = readstata13::read.dta13('.data/Engel.dta') ## get data
x.income = engel$income
y.food = engel$food

## Defining the kernel functions ##
ker <- function(v,kern){
  if (kern=="Ep") {
    K<- 0.75*(1-v^2)*(abs(v)<=1) 
  }
  if (kern=="Ga"){
    K<- 1/(sqrt(2*pi))*exp(-0.5*v^2)
  }
  return(K)
}

## Sum of kernels ##

sumkern <- function(x,x0,h,kern){
  S<- 0
  for (i in seq(1,length(x),by=1)){
    arg<-(x[i]-x0)/h
    S<-S+ker(arg,kern)
  }
  return(S)
}

## Numerator of NW estimator ##

sumkerY <- function(X,Y,x,h,kern){
  Sy<-0
  for (i in seq(1,length(X),by=1)){
    arg<-(X[i]-x)/h
    Sy<-Sy+(ker(arg,kern)*Y[i])
  }
  return(Sy)
}

## Nadaraya-Watson or a Local Constant Estimator ##
NW<- function(X,Y,x,h,kern){
  g=sumkerY(X,Y,x,h,kern)/sumkern(X,x,h,kern)	
  return(g)	
}

##### Compute bw ####

## Local Constant: LS-CS bandwidths using epanechnikov order-2
h_lc = np::npregbw(x.income,y.food,regtype="lc", bwmethod="cv.ls",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_lc = h_lc$bw

## Local Linear: LS-CS bandwidths using epanechnikov order-2
h_ll = np::npregbw(x.income,y.food,regtype="ll", bwmethod="cv.ls",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_ll =h_ll$bw

## chosen grid of points (x)
#plot(engel)
x0 = seq(200, 2000, by = 0.5)

## obtain the Local Constant estimator
g_lc <- rep(1,length(x0))
for (i in 1:length(x0)){
  g_lc[i] <- NW(x.income, y.food, x0[i], h_lc, kern = "Ep")
}
lc_res <- data.frame(x0 = x0, g_fit = g_lc) ## output

#### Plot ############
library(ggplot2)
dev.new()
ggplot2::ggplot() +
  geom_point(data = engel, aes(food,income), alpha = 0.5) +
  geom_line(data = lc_res, aes(x0, g_fit), colour = "blue") +
  xlab("Annual Household Income (in Belgian francs)") +
  ylab("Annual food expenditure") +
  ggtitle("Mean Regression Function by Local Constant")

## Local Linear Estimator ##
LL<- function(X,Y,x,h,kern){
  nn=length(X)
  e=rep(1,nn)
  X.minus.x=X-e*x
  Z=cbind(e,X.minus.x) ## constructs a n x 2 matrix  of Z.
  
  ##computing inverse matrix
  library(MASS)
  denom =matrix(0,2,2)
  for (i in seq(1,nn,by=1)){
    arg=(X[i]-x)/h
    denom=denom+ker(arg,kern)*Z[i,]%*% t(Z[i,])  ## %*% is for matrix multiplication
  } ## different from Z[i, ]*Z[i, ] component by component
  denom.inv=ginv(denom) ## to calculate the generalize inverse of matrix
  ##
  
  ##computing 2 x 1 vector
  numer=matrix(0,2,1)
  for (i in seq(1,nn,by=1)){
    arg=(X[i]-x)/h
    numer=numer+ker(arg,kern)*Z[i,]*Y[i] #to test: ker(arg,kern)*Z[1,]*Y[1]
  }
  ##
  ghat=denom.inv%*% numer
  g=ghat[1] #g(x) is the first term
  return(g)	
}

## obtain the Local Constant estimator
g_ll <- rep(1,length(x0))
for (i in 1:length(x0)){
  g_ll[i] <- LL(x.income, y.food, x0[i], h_ll, kern = "Ep")
}
ll_res <- data.frame(x0 = x0, g_fit = g_ll) ## output

#### Plot ############
library(ggplot2)
dev.new()
ggplot2::ggplot() +
  geom_point(data = engel, aes(food,income), alpha = 0.5) +
  geom_line(data = lc_res, aes(x0, g_fit), colour = "blue") +
  geom_line(data = ll_res, aes(x0, g_fit), colour = "orange")+
  xlab("Annual Household Income (in Belgian francs)") +
  ylab("Annual food expenditure") +
  ggtitle("Mean Regression Function by Local Constant")