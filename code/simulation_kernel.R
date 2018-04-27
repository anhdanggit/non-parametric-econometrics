
#----------------------------------------------------------------#
######## NONPARAMETRIC - EXAM 1: MC Simulation ###################
####### Mai-Anh Dang | Student ID: 21608631 
#
# The code includes: 
# steps of process of MC Simulation for kernel est.
#
#------------------------------------------------------------------#

library(dplyr)
library(ggplot2)
library(np)
library(magrittr)
library(readxl)
library(np)

## Define Kernel function

ker <- function(u){
  K<-1/(sqrt(2*pi))*exp(-0.5*u^2) ## [Gaussian]
  return(K)
}

## sum of kernel
sumkern <- function(x,x0,h){
  S <- 0
  for (i in seq(1,length(x),by=1)){
    arg<-(x[i]-x0)/h
    S<-S+ker(arg)
  }
  return(S)
}

#### (1) Monte-Carlo for X ###################

####### (1.1) Chosen grids of points #####################

x0 = seq(-5, 15, by = 0.01) ## chosen grids of points

## to store the summary of sim f.X
sum.fX = data.frame(n_100 = rep(0, length(x0)), n_1000 = rep(0, length(x0))) 
f.X = matrix(0, nrow = length(x0), ncol = MC) ## to store the value of each sim 
## to store the summary of sim f.X
sum.fX = data.frame(n_100 = rep(0, length(x0)), n_1000 = rep(0, length(x0))) 


######## (1.2) MC Stimulation ##############

#### START ! ###################

MC = 100
n = 100
h0 = rep(0, MC) # n = 100
h1 = rep(0, MC) # n = 1000

for (sim in (1:MC)){ ## loops of MC
  
  # random draw X
  set.seed(10*sim)
  X = rnorm(n, mean = 4, sd = 3) 
  
  # compute the h_1 (LS-CV using epanechnikov kernel)
  h_temp = npudensbw(X,bwmethod="cv.ls",bwtype="fixed", ckertype="gaussian", ckerorder=2)
  h = h_temp$bw
  #report h
  if (n == 100){
    h0[sim] = h ## store for n = 100
  } else {
    h1[sim] = h ## store for n = 1000
  }
  # estimate density function of X
  for (i in (1:length(x0)) ){ ## across the grid
    f.X[i, sim] = (1/(n*h))*sumkern(X, x0[i], h)
  }
  
  # summarize the results
  for (i in (1:length(x0))){
    if (n == 100){
    sum.fX$n_100[i] = mean(f.X[i,])
    }
    else{
      sum.fX$n_1000[i] = mean(f.X[i,])
    }
  }
}

##### (1.3) Plot ##############################

den_norm1d <- function(x, mu, sd){
  term1 = ((x-mu)^2) / (2*(sd^2))
  den = exp(-term1) / (sd*sqrt(2*pi))
  return(den)
}

mu = 4 #given
sd = 3
e = rep(1, length(x0))
den_normx <- den_norm1d(x0, e*mu, e*sd)

#### true density function

dev.new()
matplot(x0,cbind(sum.fX$n_100,sum.fX$n_1000, den_normx),type="l",col=c("red","blue", "black"),lty=c(1,1,2),
        main = expression(paste("Monte-Carlo Simulation: Density Estimation for X~N(4,9)")),
        ylab = expression(paste("Density")),
        xlab = expression(paste("X"))) 
legend("topright", legend = c("n=100", "n=1000", "true"), col=c("red", "blue", "black"), lty = c(2,1), lwd = 1)

#### report the results
res.x <- data.frame(x0 = x0, fx.n_100 = sum.fX$n_100, fx.n_1000 = sum.fX$n_1000, 
                    fx.true = den_normx)
write.csv(res.x, "result-X.csv")

#### (2) Monte-Carlo for Y ###################

## Define Sum of Multi. Kernel function

sum_mvker <- function(Y, y0_1, y0_2, h1, h2){ # Y is nx2 matrix, n 2-dim vector
  S <- 0
  y1 <- Y[,1]
  y2 <- Y[ ,2]
  for (i in (1:length(y1))){
    u1 <- (y1[i] - y0_1)/h1
    u2 <- (y2[i] - y0_2)/h2
    S <- S + ker(u1)*ker(u2)
  }
  return(S)
}


## (2.1) Chosen grids of points #####

y0_1 = seq(-5, 15, by = 0.5) ## grid in y0_1
y0_2 = seq(-5, 15, by = 0.5) ## grid in y0_2

Y0 = matrix(0, ncol = 2, nrow = length(y0_1)*length(y0_2)) # matrix 2 x points
Y0[,1] = rep(y0_1, length(y0_2))

y_02_vec = rep(y0_2[1], length(y0_1))
for (i in (2:length(y0_2))){
  y_02_vec = c(y_02_vec, rep(y0_2[i], length(y0_1)))
}

Y0[,2] = y_02_vec

# More convenient way: Y0 = expand.grid(y0_1,y0_2)

## (2.2) Simulate MC #############

# to store the values of each sim 
f.Y = matrix(0, nrow = length(y0_1)*length(y0_2), ncol = MC) 
## to store the summary of sim f.Y (joint Y1, and Y2)
sum.fY.100 = matrix(0, nrow = length(y0_1)*length(y0_2), ncol = 1) 
sum.fY.1000 = matrix(0, nrow = length(y0_1)*length(y0_2), ncol = 1) 
## to store z to plot
z = matrix(0, nrow = length(y0_1), ncol = 2)


mu = c(4,2)
sigma = matrix(c(9,1,1,4), 2, 2)

##### START! #################

MC = 100
n = 1000
sim = 1

for (sim in (1:MC)){ ## loops of MC
  
  # random draw Y
  set.seed(10*sim)
  Y = MASS::mvrnorm(n, mu, sigma)
  # random draw X
  set.seed(10*sim)
  X = rnorm(n, mean = 4, sd = 3) 
  
  # compute the h_1 (LS-CV using gaussian kernel)
  h_temp = npudensbw(X,bwmethod="cv.ls",bwtype="fixed", ckertype="gaussian", ckerorder=2)
  h = h_temp$bw
  
  # compute the h1 = h2 = h (LS-CV using epanechnikov kernel)
  h1 = h2 = h ## h0 for n = 100, h1 for n = 1000, from the previous part
  
  # calculate the join density
  for (i in (1:(length(y0_1)*length(y0_2)))){
    f.Y[i, sim] = (1/(n*h1*h2))*sum_mvker(Y, Y0[i,1], Y0[i, 2], h1, h2)
  } ## work!
  
  # summarize the results of sim
  for(i in (1:(length(y0_1)*length(y0_2)))){
    if (n == 100){
      sum.fY.100[i] = mean(f.Y[i,])
    }
    else{
      sum.fY.1000[i]= mean(f.Y[i,])
    }
  }
}

##### (2.3) Plotting ###############

out_100 <- data.frame(y0_1 = Y0 [,1], y0_2 = Y0[,2], fY.100 = sum.fY.100) ## df to plot
out_1000 <- data.frame(y0_1 = Y0 [,1], y0_2 = Y0[,2], fY.1000 = sum.fY.1000)

## as the stimulation take times, we can store the results in the csv
write.csv(out_100, file = "./results/mc_100_y.csv")
write.csv(out_1000, file = "./results/mc_1000_y.csv")

##### The true density function

den_norm2d<-function(x,y,mu1,mu2,sigma1,sigma2,rho){
  xoy = ((x-mu1)^2/sigma1^2 - 2*rho * (x-mu1)/sigma1 * (y-mu2)/sigma2 + (y-mu2)^2/sigma2^2)/(2 * (1 - rho^2))
  density = exp(-xoy)/(2 * pi *sigma1*sigma2*sqrt(1 - rho^2))
  density
}

e <- rep(1, length(Y0[,1]))
mu1 = 4 ## given
mu2 = 1
sigma1 = 3
sigma2 = 2
rho = 1/6

# true values of density
den_true <- den_norm2d(Y0[,1], Y0[,2], e*mu1, e*mu2, e*sigma1, e*sigma2, e*rho)
out_true <- data.frame(y0_1 = Y0 [,1], y0_2 = Y0[,2], fY.true = den_true) ## df to plot

### plot


dev.new()
plot.new()
scatterplot3d(out_100$y0_1, out_100$y0_2, out_100$fY.100,color="red", type = "l", pch=21,
              main="Joint Density of Y(Y1, Y2)", 
              angle = 50, box = FALSE, 
              zlab = expression(paste("Density")),
              xlab = expression(paste("Y1")),
              ylab = expression(paste("Y2")))
plot$points3d(out_true$y0_1, out_true$y0_2, out_true$fY.true, 
              col = "green", type = "l", lty = 2)
plot$points3d(out_1000$y0_1, out_1000$y0_2, out_1000$fY.1000,
              col="blue", type = "l", lty = 1)
legend("right", legend = c("n=100", "Truth", "n=1000"),
       col = c("red", "green", "blue"), pch = 16, inset = 0.1)

#### interactive plot

rgl::plot3d(out_true$y0_1, out_true$y0_2, out_true$fY.true, col = "green", type = "l",
            xlab = "Y1", ylab = "Y2", zlab = "Density")
rgl::plot3d(out_100$y0_1, out_100$y0_2, out_100$fY.100, col = "red", type = "l", add = TRUE)
rgl::plot3d(out_1000$y0_1, out_1000$y0_2, out_1000$fY.1000, col = "blue", type = "l", add = TRUE)
snapshot3d(file=paste0("3d-monto.png"))

#browseURL(paste("file://", 
#writeWebGL(dir=file.path("C:/Users/...", "3dplot-mc"), width=1400), sep=""))