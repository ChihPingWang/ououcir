rm(list=ls())
#install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
#oubmbmintegrand<-function(s, alpha.y=alpha.y){
#    alpha.y*exp(alpha.y*s)*rnorm(n=1,mean=0,sd=s)
#    }

n<-20
numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
tree<-reorder(tree,"postorder")
tree$edge
plot(tree)
nodelabels()
tiplabels()

root<-1
true.alpha.y<-1
true.alpha.tau<-0.5
true.sigma.sq.x<-2
true.sigma.tau<- 1
true.theta.tau<-3
n_s = 10
n_t = 10
true.b0 <- 0
true.b1 <- 1
true.b2 <- 0.2




##################NOT YET

model<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  alpha.tau<-model.params[2]
  sigma.sq.x<-model.params[3]
  sigma.tau<-model.params[4]
  theta.tau<-model.params[5]
  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  x1nodestates<-array(0,c(2*n-1))
  x1nodestates[n+1]<-root
  x2nodestates<-array(0,c(2*n-1))
  x2nodestates[n+1]<-root
  optimnodestates<-array(0,c(2*n-1))
  optimnodestates[n+1]<-root
  sigmasqnodestates<-array(0,c(2*n-1))
  sigmasqnodestates[n+1]<-1
  ynodestates<-array(0,c(2*n-1))
  ynodestates[n+1]<-root
  
  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length
  for(index in N:1){
    
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1nodestates[anc[index]],sd= sqrt(sigma.sq.x*treelength[index]))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2nodestates[anc[index]],sd= sqrt(sigma.sq.x*treelength[index]))
    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]
    
    c<- sigma.tau^2*(1-exp(-alpha.tau*treelength[index]))/(4*alpha.tau)
    k<- (4*theta.tau*alpha.tau)/sigma.tau^2
    lambda<-4*alpha.tau*exp(-alpha.tau*treelength[index])/(sigma.tau^2*(1-exp(-alpha.tau*treelength[index])))*sigmasqnodestates[anc[index]]
    tmp = rchisq(n=1, df=k, ncp = lambda)
    sig_u <- c*tmp
    sigmasqnodestates[des[index]]<-sig_u
    
    sigma.sq.theta<- b1^2*sigma.sq.x + b2^2*sigma.sq.x
    #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)
    
    INT1var<-treelength[index]*exp(2*alpha.y*treelength[index]) - 2*(exp(2*alpha.y*treelength[index])-exp(alpha.y*treelength[index])) + alpha.y/2*(exp(2*alpha.y*treelength[index])-1)
    INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean=0,sd=sqrt(abs(INT1var)))
    
    a <- rnorm(n=1, mean=0, sd=sqrt(theta.tau^2*(exp(2*alpha.y*treelength[index])-1)/(2*alpha.y)))
    b <- rnorm(n=1, mean=0, sd=sqrt(((sigmasqnodestates[des[index]]-theta.tau)^2/(2*(alpha.y-alpha.tau)))*(exp(2*(alpha.y-alpha.tau)*treelength[index])-1)))
    outer.int.sum=0
    for(outer.index in 1:n_t){
      inner.int.sum = 0
      for(inner.index in 1:n_s){
        c<- sigma.tau^2*(1-exp(-alpha.tau*(inner.index/n_s)))/(4*alpha.tau)
        k<- (4*theta.tau*alpha.tau)/sigma.tau^2
        lambda<- 4*alpha.tau*exp(-alpha.tau*(inner.index/n_s))/(sigma.tau^2*(1-exp(-alpha.tau*(inner.index/n_s))))*sigmasqnodestates[des[index]]
        tmp = rchisq(n=1, df=k, ncp = lambda)
        sig_u <- c*tmp
        inner.int.sum  <-  inner.int.sum + exp(alpha.tau*(inner.index/n_s))*rnorm(n=1,mean=0, sd=sqrt(1/n_s))*sqrt(sig_u)
      }
      outer.int.sum <- outer.int.sum + exp(-alpha.y*(outer.index/n_t))*inner.int.sum*rnorm(n=1,mean=0, sd=sqrt(1/n_t))
    }
    outer.int.sum
    c <- sigma.tau*outer.int.sum
    
    INTsigdWy <- exp(-alpha.y*treelength[index])*(a + b + c)
    
    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INTsigdWy
  }
  
  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  #return(c(mean(simtrait),sd(simtrait)))
  }

print(model(model.params=c(5,3,0.5,6,2),reg.params=c(2,1,3),root=root,tree=tree))

model.params=c(5,3,0.5,6,2);reg.params=c(2,1,3)
