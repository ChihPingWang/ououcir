rm(list=ls())
#install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
n<-20
numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
tree<-reorder(tree,"postorder")
tree$edge
plot(tree)
nodelabels()
tiplabels()

root<-0
true.alpha.y<-1
true.sigma.sq.y<-5
true.sigma.sq.x<-2
true.b0 <- 0
true.b1 <- 1
true.b2 <- 0.2
x1nodestates<-array(0,c(2*n-1))
x1nodestates[n+1]<-root
x2nodestates<-array(0,c(2*n-1))
x2nodestates[n+1]<-root
optimnodestates<-array(0,c(2*n-1))
optimnodestates[n+1]<-root
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
  ymean<- optimnodestates[des[index]] * (1-exp(-true.alpha.y*treelength[index])) + ynodestates[anc[index]]*exp(-true.alpha.y*treelength[index])
  sigma.sq.theta<- true.b1^2*true.sigma.sq.x + true.b2^2*true.sigma.sq.x
  yvar<-  abs((sigma.sq.theta + true.sigma.sq.y)/(2*true.alpha.y)*(1-exp(-2*true.alpha.y*treelength[index])) + sigma.sq.theta*treelength[index]*(1- 2*(1-exp(-true.alpha.y*treelength[index]))/(true.alpha.y*treelength[index])))
  ynodestates[des[index]]<-rnorm(n=1,mean=ymean,sd= sqrt(yvar))
  }

ytipstates<-ynodestates[1:n]
print(ytipstates)

ytruetrait<-ytipstates
print(ytruetrait)
summarytrait=c(mean(ytruetrait),sd(ytruetrait))

model<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  sigma.sq.y<-model.params[2]
  sigma.sq.x<-model.params[3]
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
    ymean<- optimnodestates[des[index]] * (1-exp(-alpha.y*treelength[index])) + ynodestates[anc[index]]*exp(-alpha.y*treelength[index])
    sigma.sq.theta<- b1^2*sigma.sq.x + b2^2*sigma.sq.x
    yvar<-  abs((sigma.sq.theta + sigma.sq.y)/(2*alpha.y)*(1-exp(-2*alpha.y*treelength[index])) + sigma.sq.theta*treelength[index]*(1- 2*(1-exp(-alpha.y*treelength[index]))/(alpha.y*treelength[index])))
    ynodestates[des[index]]<-rnorm(n=1,mean=ymean,sd= sqrt(yvar))
    }

  simtrait<-ynodestates[1:n]
  return(simtrait)
  #return(c(mean(simtrait),sd(simtrait)))
  }

ABC_acc<-function(model.params,reg.params,root=root,tree=tree,summarytrait=summarytrait){
  alpha.y<-model.params[1]
  sigma.sq.y<-model.params[2]
  sigma.sq.x<-model.params[3]
  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  if(alpha.y<=0 | sigma.sq.y<=0 | sigma.sq.x<=0){return(FALSE)}
  simtrait<-model(model.params,reg.params,root=root,tree=tree)
  summarysamples<-c(mean(simtrait),sd(simtrait))
  diffmean<-abs(summarysamples[1]-summarytrait[1])
  diffsd<-abs(summarysamples[2]-summarytrait[2])
  if((diffmean<0.1) & (diffsd <0.2)){
    return(TRUE)}else{return(FALSE)}
  }

run_MCMC_ABC<-function(startvalue,iterations,root=root,tree=tree,summarytrait=summarytrait){
  chain<-array(0,dim=c(iterations+1,6))
  chain[1,]<-startvalue
  for(i in 1:iterations){
    proposal<-rnorm(n=6,mean=chain[i,],sd=rep(1,6))
    model.params.proposal<-proposal[1:3]
    reg.params.proposal<-proposal[4:6]
    if(ABC_acc(model.params.proposal,reg.params.proposal, root=root, tree=tree,summarytrait=summarytrait)){
      chain[i+1,]<-proposal
    }else{
      chain[i+1,]<-chain[i,]
      }
    }
    return(mcmc(chain))
  }

posterior<-run_MCMC_ABC(c(true.alpha.y,true.sigma.sq.y,true.sigma.sq.x,true.b0,true.b1,true.b2),50000,root=root,tree=tree,summarytrait=summarytrait)
plot(posterior)
dim(posterior)
plot(mean(posterior))
apply(posterior,2,FUN=mean)
