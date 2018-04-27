#----------------------------------------------------------------#
######## NONPARAMETRIC - EXAM 2: Quantile Reg - BG90 ##########################
####### Mai-Anh Dang | Student ID: 21608631 
#
# The code includes:
# 
# estimate log(income) = g(log(food expenditure))
# (1) By Bhattacharya and Gangopadhyay (1990) Procedure
# (2) Quantile Regression
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
library(dplyr)

######### (A) BG90 ##################

## (2.1) Obtain and Explore Data ##

library(readstata13)
engel = readstata13::read.dta13('./data/Engel.dta') ## get data
engel$income = log(engel$income)
engel$food = log(engel$food)

### trimming

hist(engel$food)
engel = engel %>%
  filter(engel$food < 7.5)

# exploratory
dev.new()
rafalib::mypar(1,2)
hist(engel$food, main = paste("Histogram of X = log(expenditure)"),
     xlab = paste("log(expenditure on food)"))
hist(engel$income, main = paste("Histogram of Y = log(income)"),
     xlab = paste("log(income)"))


y.income = engel$income
x.food = engel$food

## (2.2) Kernel function ####

ker <- function(v){
   K <- 0
  if (abs(v) <= 0.5){
    K <- 1
  }
  return (K)
}

## Sum of kernels 

sumker <- function(x,x0,h){
  S = 0
  for (i in seq(1,length(x),by=1)){
    v <-(x[i]-x0)/h
    S <- S + ker(v)
  }
  return(S)
}

## sum x.y

sum.xy <- function(y, y0, x, x0, h){
  S = 0
  for (i in seq(1, length(y), by=1)){
    if (y[i] <= y0 & (abs(x[i] - x0) <= h/2)){
    S = S + 1
    }
  }
  return(S)
}

## (2.3) Chosen grids of points ####
x0 = seq(min(x.food), max(x.food), by = 0.005) 
y0 = seq(min(y.income), max(y.income), length.out = length(x0))

## (2.4) Calculate f_hat(x) ####
h_n = 0.2 ## arbitrary bw
n = length(x.food)

f_hat <- rep(0, length(x0)) # to store the value of f_hat

for (i in (1:length(x0))){
  f_hat[i]=(1/(n*h_n))*sumker(x.food,x0[i],h_n)}	

## (2.4) Calculate b_n(x) ####
b_n <- n*h_n*f_hat


## (2.5) Calculate effective alpha ####
alpha.effec = rep(0, length(b_n))
alpha.effec = sapply(0.5*b_n, as.integer) / b_n 

## (2.6) Calculate F.hat.y ####
y.x = expand.grid(y0, x0)
y.f = expand.grid(y0, f_hat)
y.alpha.effec = expand.grid(y0, alpha.effec)
y.x.f = cbind(y.x, y.f[,2])

F.hat.y = rep(0, length(y.x.f[,1]))

for (i in (1:length(y.x.f[,1]))){
  F.hat.y[i] = sum.xy(y.income, y.x.f[i,1], x.food, y.x.f[i,2], h_n)/(n*h_n*y.x.f[i,3])
}

y.x.f.F = cbind(y.x.f, F.hat.y)
y.x.f.F.alpha = cbind(y.x.f.F, y.alpha.effec[2])

result <- data.frame(y = y.x.f.F.alpha[1], x = y.x.f.F.alpha[2], F_hat = y.x.f.F.alpha[4],
                     alpha_effective = y.x.f.F.alpha[5])
colnames(result) <- c('y', 'x', 'F_hat', 'alpha_effective')

## (2.7) Take F.hat.y >= alpha_effective #####
res <- result %>%
  filter(F_hat >= alpha_effective) 

## (2.8) Fitted Y = inf #####
x =  unique(res$x)
y = rep(0, length(x))
for (i in (1: length(x))){
  y[i] = min(res$y[res$x == x[i]])
}

BG90_rel = data.frame(x0 = x, y.bg90.fit = y)
#write.csv(BG90_rel, file = "./results/bg90_est.csv")

######### (B) NPQREG ##################
temp = npqreg(txdat=as.data.frame(engel$food),
              tydat=as.data.frame(engel$income), exdat=as.data.frame(x0), tau=0.5)
fit = temp$quantile
npqreg_rel = data.frame(x0 = x0, y.npqreg.fit = fit)
#write.csv(npqreg_rel, file = "./results/npreg_est.csv")

####### (C) PLOT #######################
library(ggplot2)
library(ggthemes)
dev.new()
ggplot()+
  geom_point(data = engel, aes(x = engel$food, y = engel$income), alpha = 0.5)+
  geom_line(data = BG90_rel, aes(x = x, y = y, col = 'BG90'), size = 1)+
  geom_line(data = npqreg_rel, aes(x = x0, y = fit, col = 'npqreg'), size = 1,
            linetype = "dashed")+
  scale_colour_manual("", values = c("BG90" = "red", "npqreg"="blue"))+
  xlim(5.5, 7.5) +
  xlab("Log (Annual Household Expenditure in Food)") +
  ylab("Log (Annual Household Income)") +
  theme(legend.position = "bottom")+
  ggtitle("Quantile Regression (at median, 50th percentile)")+
  theme_bw()


  

