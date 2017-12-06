library(MASS)
rm(list=ls())
n <- 10
k <- 2
c <- 2
mu <- array(0,dim = c(k,1))
Sigma <- matrix(c(10,3,3,2),2,2)


dB_x <- mvrnorm(n=n, mu=mu, Sigm=Sigma) 
B_x1 <- cumsum(dB_x[,1])
B_x2 <- cumsum(dB_x[,2])
B_x <- cbind(B_x1,B_x2)
B_x <- rbind(0,B_x)
A <- exp(c*(1:n)/n)
B <- exp(-c*(1:n)/n)*rbind(0,dB_x[1:n-1,c(1,2)])
B1 <- cumsum(B[,1])
B2 <- cumsum(B[,2])
B <- cbind(B1,B2)
J_c <- A*B
dJ_c <- rbind(J_c[1,],diff(J_c))
int_J = t(J_c)%*%dJ_c
int_J


###############################
rm(list=ls())

n <- 10
c <- 1
mu <- 0
sigma <- 1
dB_x <- rnorm(n=n,mean = mu, sd =sigma)
B_x <- cumsum(dB_x)
B_x <- c(0,B_x)
A <- exp(c*(1:n)/n)
B <- exp(-c*(1:n)/n)*c(0,dB_x[-10])
B <- cumsum(B)
J_c <- A*B
dJ_c <- c(J_c[1],diff(J_c))
int_J = t(J_c)%*%(dJ_c)
int_J
