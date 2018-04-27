#----------------------------------------------------------------#
######## NONPARAMETRIC - EXAM 1: Density Estimation #######
####### Mai-Anh Dang | Student ID: 21608631 
#
# The code includes:
#
# estimate GDP Density by different bandwidth of kernel
#
# (1) CV Likelihood
# (2) Rule-of-thumbs
# (3) Least-squared Cross-validate
# 
# Data: Engel.dta
#------------------------------------------------------------------#

library(readxl)
library(np)

### (4.1) Obtain Data ####

gdp <- read_excel('./data/GDP.xlsx')
gdp_2005 <- gdp$`2005`
gdp_2016 <- gdp$`2016`

#gdp_2005 <- log(gdp$`2005`) ## log GDP
#gdp_2016 <- log(gdp$`2016`)

### Describe the data

pastecs::stat.desc(gdp_2005)

### (4.2) Compute the Bandwidth ####

## Rule-of-thumb
rot_2005 <- npudensbw(gdp_2005,bwmethod="normal-reference",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_rot_2005 <- rot_2005$bw
rot_2016 <- npudensbw(gdp_2016,bwmethod="normal-reference",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_rot_2016 <- rot_2016$bw

## Likelihod Cross-validation
ml_2005 <-npudensbw(gdp_2005,bwmethod="cv.ml",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_ml_2005 <- ml_2005$bw
ml_2016 <-npudensbw(gdp_2016,bwmethod="cv.ml",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_ml_2016 <- ml_2016$bw

## Least-Squared Cross-validation bw
cvls_2005 <-npudensbw(gdp_2005,bwmethod="cv.ls",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_cvls_2005 <- cvls_2005$bw
cvls_2016 <-npudensbw(gdp_2016,bwmethod="cv.ls",bwtype="fixed", ckertype="epanechnikov", ckerorder=2)
h_cvls_2016 <- cvls_2016$bw

## Report
bw.output<- c("GDP_year","rule-of-thumb","cv-ml","ls-cv","2005",h_rot_2005, h_ml_2005, h_cvls_2005,
              "2016",h_rot_2016, h_ml_2016, h_cvls_2016)
bw.mat<-matrix(bw.output,nrow=3,ncol=4,byrow=TRUE)  
rm(bw.output)
print(bw.mat)
bw_report <- data.frame(Year = c(2005, 2016), 
                        h_ml = c(h_ml_2005, h_ml_2016),
                        h_rot = c(h_rot_2005, h_rot_2016),
                        h_cvls = c(h_cvls_2005, h_cvls_2016)
)
print(xtable::xtable(bw_report, type = "latex", file = "filename2.tex"))

## (4.3) Estimate density function ####

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

### choosing a grid of points

z <- seq(0,1500,by=0.5) 
#z <- seq(-3, 10, 0.005) ### log GDP 

f.2005<-matrix(0,nrow=length(z),ncol=3)
f.2016<-matrix(0,nrow=length(z),ncol=3)

n <- length(gdp_2005)

for (i in seq(1,length(z),by=1)){
  f.2005[i,1]=(1/(n*h_ml_2005))*sumkern(gdp_2005,z[i],h_ml_2005,"Ep")	
  f.2005[i,2]=(1/(n*h_rot_2005))*sumkern(gdp_2005,z[i],h_rot_2005,"Ep")	
  f.2005[i,3]=(1/(n*h_cvls_2005))*sumkern(gdp_2005,z[i],h_cvls_2005,"Ep")	
  
  f.2016[i,1]=(1/(n*h_ml_2016))*sumkern(gdp_2016,z[i],h_ml_2016,"Ep")	
  f.2016[i,2]=(1/(n*h_rot_2016))*sumkern(gdp_2005,z[i],h_rot_2016,"Ep")	
  f.2016[i,3]=(1/(n*h_cvls_2016))*sumkern(gdp_2005,z[i],h_cvls_2016,"Ep")		
}

###### (4.4) Plot ######

dev.new()
matplot(z,cbind(f.2005[,1],f.2016[,1]),type="l",col=c("red","blue"),lty=c(2,1),
        main = expression(paste("Density Estimation by Likelihood Cross-Validation")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log GDP (in billions of dollars)"))) 
legend("topright", legend = c("2005", "2016"), col=c("red", "blue"), lty = c(2,1), lwd = 1)

dev.new()
matplot(z,cbind(f.2005[,2],f.2016[,2]),type="l",col=c("red","blue"),lty=c(2,1),
        main = expression(paste("Density Estimation by Rule-of-Thumbs")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log GDP (in billions of dollars)"))) 
legend("topright", legend = c("2005", "2016"), col=c("red", "blue"), lty = c(1,1), lwd = 1)

dev.new()
matplot(z,cbind(f.2005[,3],f.2016[,3]),type="l",col=c("red","blue"),lty=c(2,1),
        main = expression(paste("Density Estimation by Least-Squared Cross-Validation Bandwidth")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("Log GDP (in billions of dollars)"))) 
legend("topright", legend = c("2005", "2016"), col=c("red", "blue"), lty = c(2,1), lwd = 1)  

