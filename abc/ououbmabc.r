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

root<-0
true.alpha.y<-1
true.alpha.x<-0.2
true.optim.x<-0.5
true.sigma.sq.x<-2
true.tau<- 1

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
sigmanodestates[n+1]<-root
ynodestates<-array(0,c(2*n-1))
ynodestates[n+1]<-root

N<-dim(tree$edge)[1]
anc<-tree$edge[,1]
des<-tree$edge[,2]
treelength<-tree$edge.length

int3integrand<- function(s,alpha.x=alpha.x,sigma.sq.theta=sigma.sq.theta,alpha.y=alpha.y)
  int3stochV<-(exp(2* alpha.x) -1 ) / (2*alpha.x)
  f<-sqrt(sigma.sq.theta)*alpha.y*exp((alpha.y-alpha.x)*t)*1/sqrt(2*pi*int3stochV)*exp(-0.5/int3stochV*s^2)
  return(f)
  }
for(index in N:1){
  x1ou.mean<-x1nodestates[anc[index]]*exp(-true.alpha.x*treelength[index]) + true.optim.x*(1-exp(-true.alpha.x*treelength[index]))
  x1ou.sd<-sqrt((true.sigma.sq.x/(2*true.alpha.x))*(1-exp(-2*true.alpha.x*treelength[index])))
  x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd=x1ou.sd)

  x2ou.mean<-x2nodestates[anc[index]]*exp(-true.alpha.x*treelength[index]) + true.optim.x*(1-exp(-true.alpha.x*treelength[index]))
  x2ou.sd<-sqrt((true.sigma.sq.x/(2*true.alpha.x))*(1-exp(-2*true.alpha.x*treelength[index])))
  x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd=x2ou.sd)

  optimnodestates[des[index]]<- true.b0 + true.b1*x1nodestates[des[index]] + true.b2*x2nodestates[des[index]]

  sigmasqnodestates[des[index]]<-rnorm(n=1,mean=sigmasqnodestates[anc[index]],sd=sqrt(true.tau*treelength[index]))
  sigma.sq.theta<- true.b1^2*true.sigma.sq.x + true.b2^2*true.sigma.sq.x
  #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)

  A1<-  (true.alpha.y*optimnodestates[des[index]]/ (true.alpha.y-true.alpha.x)) *(exp((true.alpha.y-true.alpha.x)*treelength[index]) -1)
  theta1<-0
  A2<- theta1*(exp(true.alpha.y*treelength[index])-1) - (true.alpha.y*theta1/(true.alpha.y-true.alpha.x))*(exp((true.alpha.y-true.alpha.x)*treelength[index]) -1)
  # WRONG A3<- integrate(int3intgrand,lower=0,upper=treelength[index],alpha.x=true.alpha.x,sigma.sq.theta=true.sigma.sq.theta,alpha.y=true.alpha.y)
  A3<-1/(sigma.sq.theta*true.alpha.y^2)*rnorm(1,mean=0,sd=sqrt( (sigma.sq.theta*true.alpha.y^2)^3 *treelength[index]^3 /3    ))
  INTtime<- A1+A2+A3
  INT1<-exp(-true.alpha.y*treelength[index])*INTtime

  fexpr<-expression(exp(true.alpha.y*t)*w)
  res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
  INT2<-exp(-true.alpha.y*treelength[index])*median(res$X)

  ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
  }

ytipstates<-ynodestates[1:n]
print(ytipstates)

ytruetrait<-ytipstates
print(ytruetrait)
summarytrait=c(mean(ytruetrait),sd(ytruetrait))


##################NOT YET

model<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  alpha.x<-model.params[2]
  optima.x<-model.params[3]
  sigm.sq.x<-model.params[4]
  tau<-model.params[5]

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
  sigmanodestates[n+1]<-root
  ynodestates<-array(0,c(2*n-1))
  ynodestates[n+1]<-root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length

  for(index in N:1){
    x1ou.mean<-x1nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + optim.x*(1-exp(-alpha.x*treelength[index]))
    x1ou.sd<-sqrt((sigma.sq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd=x1ou.sd)

    x2ou.mean<-x2nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + optim.x*(1-exp(-alpha.x*treelength[index]))
    x2ou.sd<-sqrt((sigma.sq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd=x2ou.sd)

    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]

    sigmasqnodestates[des[index]]<-rnorm(n=1,mean=sigmasqnodestates[anc[index]],sd=sqrt(tau*treelength[index]))
    sigma.sq.theta<- b1^2*sigma.sq.x + b2^2*sigma.sq.x
    #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)

    A1<-  (alpha.y*optimnodestates[des[index]]/ (alpha.y-alpha.x)) *(exp((alpha.y-alpha.x)*treelength[index]) -1)
    theta1<-0
    A2<- theta1*(exp(alpha.y*treelength[index])-1) - (alpha.y*theta1/(alpha.y-alpha.x))*(exp((alpha.y-alpha.x)*treelength[index]) -1)
    A3<-1/(sigma.sq.theta*alpha.y^2)*rnorm(1,mean=0,sd=sqrt((sigma.sq.theta*alpha.y^2)^3 *tree[index]^3 /3))
    INTtime<- A1+A2+A3
    INT1<-exp(-alpha.y*treelength[index])*INTtime

    fexpr<-expression(alpha.y*exp(alpha.y*t)*w)
    res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
    INT2<-median(res$X)

    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
    }

  simtrait<-ynodestates[1:n]
  return(c(mean(simtrait),sd(simtrait)))
  }

ABC_acc<-function(model.params,reg.params,root=root,tree=tree,summarytrait=summarytrait){
  alpha.y<-model.params[1]
  alpha.x<-model.params[2]
  optima.x<-model.params[3]
  sigm.sq.x<-model.params[4]
  tau<-model.params[5]

  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  if(alpha.y<=0 | alpha.x<=0 |sigma.sq.x<=0 | tau<=0){return(FALSE)}
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
    proposal<-rnorm(n=8,mean=chain[i,],sd=rep(1,8))
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

posterior<-run_MCMC_ABC(c(true.alpha.y,true.alpha.x,true.optim.x,true.sigma.sq.x,true.tau,true.b0,true.b1,true.b2),500000,root=root,tree=tree,summarytrait=summarytrait)
plot(posterior)
plot(mean(posterior))
