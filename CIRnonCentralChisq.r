rm(list = ls())
sig_tau=2
alp_y=1
alp_tau=5
theta_tau=4
u=1
sig_0=3
n_s = 1000
n_t = 1000
t=1
# c<- sig_tau^2*(1-exp(-alp_tau*u))/(4*alp_tau)
# k<- (4*theta_tau*alp_tau)/sig_tau^2
# lambda<-4*alp_tau*exp(-alp_tau*u)/(sig_tau^2*(1-exp(-alp_tau*u)))*y_0 
# tmp = rchisq(n=1, df=k, ncp = lambda) 
# sig_u <- c*tmp
# sig_u
# sqrt(sig_u)
###(C)###
outer.int.sum=0
for(outer.index in 1:n_t){
    inner.int.sum = 0
    for(inner.index in 1:n_s){
      c<- sig_tau^2*(1-exp(-alp_tau*(inner.index/n_s)))/(4*alp_tau)
      k<- (4*theta_tau*alp_tau)/sig_tau^2
      lambda<- 4*alp_tau*exp(-alp_tau*(inner.index/n_s))/(sig_tau^2*(1-exp(-alp_tau*(inner.index/n_s))))*sig_0
      tmp = rchisq(n=1, df=k, ncp = lambda)
      sig_u <- c*tmp
      inner.int.sum  <-  inner.int.sum + exp(alp_tau*(inner.index/n_s))*rnorm(n=1,mean=0, sd=sqrt(1/n_s))*sqrt(sig_u)
    }
  outer.int.sum <- outer.int.sum + exp(-alp_y*(outer.index/n_t))*inner.int.sum*rnorm(n=1,mean=0, sd=sqrt(1/n_t))
}
outer.int.sum
c <- sig_tau*outer.int.sum

###(B)###
b <- rnorm(n=1, mean=0, sd=sqrt(((sig_0-theta_tau)^2/(2*(alp_y-alp_tau)))*(exp(2*(alp_y-alp_tau)*t)-1)))

###(A)###
a <- rnorm(n=1, mean=0, sd=sqrt(theta_tau^2*(exp(2*alp_y*t)-1)/(2*alp_y)))
a + b + c
