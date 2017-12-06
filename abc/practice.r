install.packages("EasyABC")
rm(list=ls())
library(EasyABC)
#my_prior=list(c("unif",0,1),c("normal",1,2))
#my_prior=list(list(c("runif",1,0,1), c("dunif",0,1)))
toy_model<-function(x){
  c(x[1]+x[2]+rnorm(1,0,0.1),x[1]*x[2]+rnorm(1,0,0.1))
}
toy_prior = list(c("unif",0,1),c("normal",1,2))
sum_stat_obs=c(1.5,0.5)
set.seed(1)
n=10
p=0.5
ABC_rej<-ABC_rejection(model=toy_model,prior=toy_prior,nb_simul=n,summary_stat_target=sum_stat_obs,tol=p)

?ABC_rejection
#https://stats.stackexchange.com/questions/277499/simple-linear-regression-using-approximate-bayesian-computation-abc/277512