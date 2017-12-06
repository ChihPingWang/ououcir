rm(list=ls())
#install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
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
  
  x1nodestates[des[index]]<-rnorm(n=1,mean=x1nodestates[anc[index]],sd= sqrt(true.sigma.sq.x*treelength[index]))
  x2nodestates[des[index]]<-rnorm(n=1,mean=x2nodestates[anc[index]],sd= sqrt(true.sigma.sq.x*treelength[index]))
  optimnodestates[des[index]]<- true.b0 + true.b1*x1nodestates[des[index]] + true.b2*x2nodestates[des[index]]

  c<- true.sigma.tau^2*(1-exp(-true.alpha.tau*treelength[index]))/(4*true.alpha.tau)
  k<- (4*true.theta.tau*true.alpha.tau)/true.sigma.tau^2
  lambda<-4*true.alpha.tau*exp(-true.alpha.tau*treelength[index])/(true.sigma.tau^2*(1-exp(-true.alpha.tau*treelength[index])))*sigmasqnodestates[anc[index]]
  tmp = rchisq(n=1, df=k, ncp = lambda)
  sig_u <- c*tmp
  sigmasqnodestates[des[index]]<-sig_u

  sigma.sq.theta<- true.b1^2*true.sigma.sq.x + true.b2^2*true.sigma.sq.x
  #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)

  INT1var<-treelength[index]*exp(2*true.alpha.y*treelength[index]) - 2*(exp(2*true.alpha.y*treelength[index])-exp(true.alpha.y*treelength[index])) + true.alpha.y/2*(exp(2*true.alpha.y*treelength[index])-1)
  INT1<-exp(-true.alpha.y*treelength[index])* rnorm(n=1,mean=0,sd=sqrt(abs(INT1var)))

  a <- rnorm(n=1, mean=0, sd=sqrt(true.theta.tau^2*(exp(2*true.alpha.y*treelength[index])-1)/(2*true.alpha.y)))
  b <- rnorm(n=1, mean=0, sd=sqrt(((sigmasqnodestates[des[index]]-true.theta.tau)^2/(2*(true.alpha.y-true.alpha.tau)))*(exp(2*(true.alpha.y-true.alpha.tau)*treelength[index])-1)))
  outer.int.sum=0
  for(outer.index in 1:n_t){
      inner.int.sum = 0
      for(inner.index in 1:n_s){
        c<- true.sigma.tau^2*(1-exp(-true.alpha.tau*(inner.index/n_s)))/(4*true.alpha.tau)
        k<- (4*true.theta.tau*true.alpha.tau)/true.sigma.tau^2
        lambda<- 4*true.alpha.tau*exp(-true.alpha.tau*(inner.index/n_s))/(true.sigma.tau^2*(1-exp(-true.alpha.tau*(inner.index/n_s))))*sigmasqnodestates[des[index]]
        tmp = rchisq(n=1, df=k, ncp = lambda)
        sig_u <- c*tmp
        inner.int.sum  <-  inner.int.sum + exp(true.alpha.tau*(inner.index/n_s))*rnorm(n=1,mean=0, sd=sqrt(1/n_s))*sqrt(sig_u)
      }
    outer.int.sum <- outer.int.sum + exp(-true.alpha.y*(outer.index/n_t))*inner.int.sum*rnorm(n=1,mean=0, sd=sqrt(1/n_t))
  }
  outer.int.sum
  c <- true.sigma.tau*outer.int.sum

  INTsigdWy <- exp(-true.alpha.y*treelength[index])*(a + b + c)

  ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INTsigdWy
  }

ytipstates<-ynodestates[1:n]
print(ytipstates)

ytruetrait<-ytipstates
print(ytruetrait)
summarytrait=c(mean(ytruetrait),sd(ytruetrait))


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

  n<-Ntip(tree)
  x1nodestates<-array(0,c(2*n-1))
  x1nodestates[n+1]<-root
  x2nodestates<-array(0,c(2*n-1))
  x2nodestates[n+1]<-root
  optimnodestates<-array(0,c(2*n-1))
  optimnodestates[n+1]<-root
  sigmasqnodestates<-array(0,c(2*n-1))
  sigmasqnodestates[n+1]<-root
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
    sigmasqnodestates[des[index]]<-rnorm(n=1,mean=sigmasqnodestates[anc[index]],sd=sqrt(sigma.tau*treelength[index]))
    sigma.sq.theta<- b1^2*sigma.sq.x + b2^2*sigma.sq.x
    #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)

    INT1var<-treelength[index]*exp(2*alpha.y*treelength[index]) - 2*(exp(2*alpha.y*treelength[index])-exp(alpha.y*treelength[index])) + alpha.y/2*(exp(2*alpha.y*treelength[index])-1)
    INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean=0,sd=sqrt(sigma.sq.theta*abs(INT1var)))
    print(sqrt(sigma.sq.theta*abs(INT1var)))### debug
    
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
    
    INTsigdWy <- a + b + c
    

    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INTsigdWy 
    }

  simtrait<-ynodestates[1:n]
  return(summarytrait=c(mean(simtrait),sd(simtrait)))
  }

ABC_acc<-function(model.params,reg.params,root=root,tree=tree,summarytrait=summarytrait){
  alpha.y<-model.params[1]
  alpha.tau<-model.params[2]
  sigma.sq.x<-model.params[3]
  sigma.tau<-model.params[4]
  theta.tau<-model.params[5]
  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  if(alpha.y<=0 | alpha.tau<=0 | sigma.sq.x<=0 | sigma.tau<=0 | theta.tau<=0){return(FALSE)}
  summarysamples<-model(model.params,reg.params,root=root,tree=tree)
  diffmean<-abs(summarysamples[1]-summarytrait[1])
  diffsd<-abs(summarysamples[2]-summarytrait[2])
  if((diffmean<0.1) & (diffsd <0.2)){
    return(TRUE)}else{return(FALSE)}
  }

run_MCMC_ABC<-function(startvalue,iterations,root=root,tree=tree,summarytrait=summarytrait){
  chain<-array(0,dim=c(iterations+1,8))
  chain[1,]<-startvalue
  for(i in 1:iterations){
    proposal<-rnorm(n=6,mean=chain[i,],sd=rep(1,8))
    model.params.proposal<-proposal[1:5]
    reg.params.proposal<-proposal[6:8]
    if(ABC_acc(model.params.proposal,reg.params.proposal, root=root, tree=tree,summarytrait=summarytrait)){
      chain[i+1,]<-proposal
    }else{
      chain[i+1,]<-chain[i,]
      }
    }
    return(mcmc(chain))
  }

posterior<-run_MCMC_ABC(c(true.alpha.y,true.alpha.tau,true.sigma.sq.x,true.sigma.tau,true.theta.tau,true.b0,true.b1,true.b2),5000,root=root,tree=tree,summarytrait=summarytrait)
plot(posterior)
plot(mean(posterior))
