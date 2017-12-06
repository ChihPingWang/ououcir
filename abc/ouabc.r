rm(list=ls())
#install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
n<-20
numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
?TreeSim
tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
tree<-reorder(tree,"postorder")
tree$edge
plot(tree)
nodelabels()
tiplabels()

root<-0
true.alpha<-1
true.optim<-5
true.sigma.sq<-2
nodestates<-array(0,c(2*n-1))
nodestates[n+1]<-root
N<-dim(tree$edge)[1]
anc<-tree$edge[,1]
des<-tree$edge[,2]
treelength<-tree$edge.length
for(index in N:1){

#  print(c(anc[index],des[index]))
  ou.mean<-nodestates[anc[index]]*exp(-true.alpha*treelength[index]) + true.optim*(1-exp(-true.alpha*treelength[index]))
  ou.sd<-sqrt((true.sigma.sq/(2*true.alpha))*(1-exp(-2*true.alpha*treelength[index])))
  nodestates[des[index]]<-rnorm(n=1,mean=ou.mean,sd=ou.sd )
  }
tipstates<-nodestates[1:n]
print(tipstates)

truetrait<-tipstates
summarytrait=c(mean(truetrait),sd(truetrait))



model<-function(model.params,root=root,tree=tree){
  alpha<-model.params[1]
  mu<-model.params[2]
  sigma.sq<-model.params[3]
  #sigmasq<-1
  n<-Ntip(tree)
  nodestates<-array(0,c(2*n-1))
  nodestates[n+1]<-root
  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length

  for(index in N:1){
    ou.mean<-nodestates[anc[index]]*exp(-alpha*treelength[index]) + mu*(1-exp(-alpha*treelength[index]))
    ou.sd<-sqrt((sigma.sq/(2*alpha))*(1-exp(-2*alpha*treelength[index])))
    nodestates[des[index]]<-rnorm(n=1,mean=ou.mean,sd=ou.sd)
    }
  simtrait<-nodestates[1:n]
  return(c(mean(simtrait),sd(simtrait)))
}


ABC_acc<-function(model.params,root=root,tree=tree,summarytrait=summarytrait){
    alpha<-model.params[1]
    mu<-model.params[2]
    sigma.sq<-model.params[3]
    if(alpha<=0  | sigma.sq<=0){return(FALSE)}
    summarysamples<-model(model.params,root=root,tree=tree)
    diffmean<-abs(summarysamples[1]-summarytrait[1])
    diffsd<-abs(summarysamples[2]-summarytrait[2])
    if( (diffmean<0.1) & (diffsd<0.2)){
      return(TRUE)}else{return(FALSE)}
      }


run_MCMC_ABC<-function(startvalue,iterations,root=root,tree=tree,summarytrait=summarytrait){
  chain<-array(0,dim=c(iterations+1,3))
  chain[1,]<-startvalue
  for(i in 1:iterations){
    proposal<-rnorm(n=3,mean=chain[i,],sd=c(1,1,1))
    if(ABC_acc(proposal, root=root,tree=tree,summarytrait=summarytrait)){
    chain[i+1,]<-proposal
    }else{
    chain[i+1,]<-chain[i,]
    }
  }
  return(mcmc(chain))
  }

posterior<-run_MCMC_ABC(c(1,5,2),50000,root=root,tree=tree,summarytrait=summarytrait)
plot(posterior)
print(mean(posterior))
