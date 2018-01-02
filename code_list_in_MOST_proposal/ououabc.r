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
true.alpha.x<-0.2
true.optim.x<-0.5
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
  x1ou.mean<-x1nodestates[anc[index]]*exp(-true.alpha.x*treelength[index]) + true.optim.x*(1-exp(-true.alpha.x*treelength[index]))
  x1ou.sd<-sqrt((true.sigma.sq.x/(2*true.alpha.x))*(1-exp(-2*true.alpha.x*treelength[index])))
  x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd= x1ou.sd)

  x2ou.mean<-x2nodestates[anc[index]]*exp(-true.alpha.x*treelength[index]) + true.optim.x*(1-exp(-true.alpha.x*treelength[index]))
  x2ou.sd<-sqrt((true.sigma.sq.x/(2*true.alpha.x))*(1-exp(-2*true.alpha.x*treelength[index])))
  x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd= x2ou.sd)

  optimnodestates[des[index]]<- true.b0 + true.b1*x1nodestates[des[index]] + true.b2*x2nodestates[des[index]]
  sigma.sq.theta<- true.b1^2*true.sigma.sq.x + true.b2^2*true.sigma.sq.x
  ymean<- ((true.alpha.y - true.alpha.x)*ynodestates[anc[index]]*exp(true.alpha.x*treelength[index]) + (true.alpha.y*exp(true.alpha.y*treelength[index]) - true.alpha.y*exp(true.alpha.x*treelength[index]))*optimnodestates[des[index]])*exp(-true.alpha.y*treelength[index] - true.alpha.x*treelength[index])/(true.alpha.y - true.alpha.x)
  #E_y_sol <- ((alp - bet)*y_0*e^(bet*ta) + (alp*e^(alp*ta) - alp*e^(bet*ta))*th_0)*e^(-alp*ta - bet*ta)/(alp - bet)

  rho_y_th <- 0
  yvar<- -1/2*((true.alpha.y^4 + true.alpha.y^3*true.alpha.x)*sigma.sq.theta^2*exp(2*true.alpha.y*treelength[index]) - 2*(true.alpha.y^4*true.alpha.x - true.alpha.y^3*true.alpha.x^2 - true.alpha.y^2*true.alpha.x^3 + true.alpha.y*true.alpha.x^4)*ynodestates[anc[index]]^2*exp(2*true.alpha.x*treelength[index]) - 2*((true.alpha.y^4*true.alpha.x + true.alpha.y^3*true.alpha.x^2)*exp(2*true.alpha.y*treelength[index]) - 2*(true.alpha.y^4*true.alpha.x + true.alpha.y^3*true.alpha.x^2)*exp(true.alpha.y*treelength[index] + true.alpha.x*treelength[index]) + (true.alpha.y^4*true.alpha.x + true.alpha.y^3*true.alpha.x^2)*exp(2*true.alpha.x*treelength[index]))*optimnodestates[des[index]]^2 - 4*((true.alpha.y^4*true.alpha.x - true.alpha.y^2*true.alpha.x^3)*exp(true.alpha.y*treelength[index] + true.alpha.x*treelength[index]) - (true.alpha.y^4*true.alpha.x - true.alpha.y^2*true.alpha.x^3)*exp(2*true.alpha.x*treelength[index]))*optimnodestates[des[index]]*ynodestates[anc[index]] - 4*(true.alpha.y^3*true.alpha.x*sigma.sq.theta^2 - (true.alpha.y^3*true.alpha.x - true.alpha.y^2*true.alpha.x^2)*sigma.sq.theta*rho_y_th*true.sigma.sq.y)*exp(true.alpha.y*treelength[index] + true.alpha.x*treelength[index]) - (2*(true.alpha.y^3*true.alpha.x - true.alpha.y*true.alpha.x^3)*sigma.sq.theta*rho_y_th*true.sigma.sq.y - (true.alpha.y^3*true.alpha.x + true.alpha.y^2*true.alpha.x^2)*sigma.sq.theta^2 - (true.alpha.y^3*true.alpha.x - true.alpha.y^2*true.alpha.x^2 - true.alpha.y*true.alpha.x^3 + true.alpha.x^4)*true.sigma.sq.y^2 + (2*(true.alpha.y^3*true.alpha.x - 2*true.alpha.y^2*true.alpha.x^2 + true.alpha.y*true.alpha.x^3)*sigma.sq.theta*rho_y_th*true.sigma.sq.y + (true.alpha.y^4 - 2*true.alpha.y^3*true.alpha.x + true.alpha.y^2*true.alpha.x^2)*sigma.sq.theta^2 + (true.alpha.y^3*true.alpha.x - true.alpha.y^2*true.alpha.x^2 - true.alpha.y*true.alpha.x^3 + true.alpha.x^4)*true.sigma.sq.y^2)*exp(2*true.alpha.y*treelength[index]))*exp(2*true.alpha.x*treelength[index]))*exp(-2*true.alpha.y*treelength[index] - 2*true.alpha.x*treelength[index])/(true.alpha.y^4*true.alpha.x - true.alpha.y^3*true.alpha.x^2 - true.alpha.y^2*true.alpha.x^3 + true.alpha.y*true.alpha.x^4) - ((true.alpha.y - true.alpha.x)*ynodestates[anc[index]]*exp(true.alpha.x*treelength[index]) + (true.alpha.y*exp(true.alpha.y*treelength[index]) - true.alpha.y*exp(true.alpha.x*treelength[index]))*optimnodestates[des[index]])^2*exp(-2*true.alpha.y*treelength[index] - 2*true.alpha.x*treelength[index])/(true.alpha.y - true.alpha.x)^2

  #var_y<- -1/2*((alp^4 + alp^3*bet)*eta^2*e^(2*alp*ta) - 2*(alp^4*bet - alp^3*bet^2 - alp^2*bet^3 + alp*bet^4)*y_0^2*e^(2*bet*ta) - 2*((alp^4*bet + alp^3*bet^2)*e^(2*alp*ta) - 2*(alp^4*bet + alp^3*bet^2)*e^(alp*ta + bet*ta) + (alp^4*bet + alp^3*bet^2)*e^(2*bet*ta))*th_0^2 - 4*((alp^4*bet - alp^2*bet^3)*e^(alp*ta + bet*ta) - (alp^4*bet - alp^2*bet^3)*e^(2*bet*ta))*th_0*y_0 - 4*(alp^3*bet*eta^2 - (alp^3*bet - alp^2*bet^2)*eta*rho_y_th*sig)*e^(alp*ta + bet*ta) - (2*(alp^3*bet - alp*bet^3)*eta*rho_y_th*sig - (alp^3*bet + alp^2*bet^2)*eta^2 - (alp^3*bet - alp^2*bet^2 - alp*bet^3 + bet^4)*sig^2 + (2*(alp^3*bet - 2*alp^2*bet^2 + alp*bet^3)*eta*rho_y_th*sig + (alp^4 - 2*alp^3*bet + alp^2*bet^2)*eta^2 + (alp^3*bet - alp^2*bet^2 - alp*bet^3 + bet^4)*sig^2)*e^(2*alp*ta))*e^(2*bet*ta))*e^(-2*alp*ta - 2*bet*ta)/(alp^4*bet - alp^3*bet^2 - alp^2*bet^3 + alp*bet^4) - ((alp - bet)*y_0*e^(bet*ta) + (alp*e^(alp*ta) - alp*e^(bet*ta))*th_0)^2*e^(-2*alp*ta - 2*bet*ta)/(alp - bet)^2

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
  alpha.x<-model.params[3]
  optim.x<-model.params[4]
  sigma.sq.x<-model.params[5]

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
    x1ou.mean<-x1nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + optim.x*(1-exp(-alpha.x*treelength[index]))
    x1ou.sd<-sqrt((sigma.sq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd= x1ou.sd)

    x2ou.mean<-x2nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + optim.x*(1-exp(-true.alpha.x*treelength[index]))
    x2ou.sd<-sqrt((sigma.sq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd= x2ou.sd)

    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]
    sigma.sq.theta<- b1^2*sigma.sq.x + b2^2*sigma.sq.x

    ymean<- ((true.alpha.y - true.alpha.x)*ynodestates[anc[index]]*exp(true.alpha.x*treelength[index]) + (true.alpha.y*exp(true.alpha.y*treelength[index]) - true.alpha.y*exp(true.alpha.x*treelength[index]))*optimnodestates[des[index]])*exp(-true.alpha.y*treelength[index] - true.alpha.x*treelength[index])/(true.alpha.y - true.alpha.x)
    #E_y_sol <- ((alp - bet)*y_0*e^(bet*ta) + (alp*e^(alp*ta) - alp*e^(bet*ta))*th_0)*e^(-alp*ta - bet*ta)/(alp - bet)

    rho_y_th <- 0
    yvar<- -1/2*((true.alpha.y^4 + true.alpha.y^3*true.alpha.x)*sigma.sq.theta^2*exp(2*true.alpha.y*treelength[index]) - 2*(true.alpha.y^4*true.alpha.x - true.alpha.y^3*true.alpha.x^2 - true.alpha.y^2*true.alpha.x^3 + true.alpha.y*true.alpha.x^4)*ynodestates[anc[index]]^2*exp(2*true.alpha.x*treelength[index]) - 2*((true.alpha.y^4*true.alpha.x + true.alpha.y^3*true.alpha.x^2)*exp(2*true.alpha.y*treelength[index]) - 2*(true.alpha.y^4*true.alpha.x + true.alpha.y^3*true.alpha.x^2)*exp(true.alpha.y*treelength[index] + true.alpha.x*treelength[index]) + (true.alpha.y^4*true.alpha.x + true.alpha.y^3*true.alpha.x^2)*exp(2*true.alpha.x*treelength[index]))*optimnodestates[des[index]]^2 - 4*((true.alpha.y^4*true.alpha.x - true.alpha.y^2*true.alpha.x^3)*exp(true.alpha.y*treelength[index] + true.alpha.x*treelength[index]) - (true.alpha.y^4*true.alpha.x - true.alpha.y^2*true.alpha.x^3)*exp(2*true.alpha.x*treelength[index]))*optimnodestates[des[index]]*ynodestates[anc[index]] - 4*(true.alpha.y^3*true.alpha.x*sigma.sq.theta^2 - (true.alpha.y^3*true.alpha.x - true.alpha.y^2*true.alpha.x^2)*sigma.sq.theta*rho_y_th*true.sigma.sq.y)*exp(true.alpha.y*treelength[index] + true.alpha.x*treelength[index]) - (2*(true.alpha.y^3*true.alpha.x - true.alpha.y*true.alpha.x^3)*sigma.sq.theta*rho_y_th*true.sigma.sq.y - (true.alpha.y^3*true.alpha.x + true.alpha.y^2*true.alpha.x^2)*sigma.sq.theta^2 - (true.alpha.y^3*true.alpha.x - true.alpha.y^2*true.alpha.x^2 - true.alpha.y*true.alpha.x^3 + true.alpha.x^4)*true.sigma.sq.y^2 + (2*(true.alpha.y^3*true.alpha.x - 2*true.alpha.y^2*true.alpha.x^2 + true.alpha.y*true.alpha.x^3)*sigma.sq.theta*rho_y_th*true.sigma.sq.y + (true.alpha.y^4 - 2*true.alpha.y^3*true.alpha.x + true.alpha.y^2*true.alpha.x^2)*sigma.sq.theta^2 + (true.alpha.y^3*true.alpha.x - true.alpha.y^2*true.alpha.x^2 - true.alpha.y*true.alpha.x^3 + true.alpha.x^4)*true.sigma.sq.y^2)*exp(2*true.alpha.y*treelength[index]))*exp(2*true.alpha.x*treelength[index]))*exp(-2*true.alpha.y*treelength[index] - 2*true.alpha.x*treelength[index])/(true.alpha.y^4*true.alpha.x - true.alpha.y^3*true.alpha.x^2 - true.alpha.y^2*true.alpha.x^3 + true.alpha.y*true.alpha.x^4) - ((true.alpha.y - true.alpha.x)*ynodestates[anc[index]]*exp(true.alpha.x*treelength[index]) + (true.alpha.y*exp(true.alpha.y*treelength[index]) - true.alpha.y*exp(true.alpha.x*treelength[index]))*optimnodestates[des[index]])^2*exp(-2*true.alpha.y*treelength[index] - 2*true.alpha.x*treelength[index])/(true.alpha.y - true.alpha.x)^2

    #var_y<- -1/2*((alp^4 + alp^3*bet)*eta^2*e^(2*alp*ta) - 2*(alp^4*bet - alp^3*bet^2 - alp^2*bet^3 + alp*bet^4)*y_0^2*e^(2*bet*ta) - 2*((alp^4*bet + alp^3*bet^2)*e^(2*alp*ta) - 2*(alp^4*bet + alp^3*bet^2)*e^(alp*ta + bet*ta) + (alp^4*bet + alp^3*bet^2)*e^(2*bet*ta))*th_0^2 - 4*((alp^4*bet - alp^2*bet^3)*e^(alp*ta + bet*ta) - (alp^4*bet - alp^2*bet^3)*e^(2*bet*ta))*th_0*y_0 - 4*(alp^3*bet*eta^2 - (alp^3*bet - alp^2*bet^2)*eta*rho_y_th*sig)*e^(alp*ta + bet*ta) - (2*(alp^3*bet - alp*bet^3)*eta*rho_y_th*sig - (alp^3*bet + alp^2*bet^2)*eta^2 - (alp^3*bet - alp^2*bet^2 - alp*bet^3 + bet^4)*sig^2 + (2*(alp^3*bet - 2*alp^2*bet^2 + alp*bet^3)*eta*rho_y_th*sig + (alp^4 - 2*alp^3*bet + alp^2*bet^2)*eta^2 + (alp^3*bet - alp^2*bet^2 - alp*bet^3 + bet^4)*sig^2)*e^(2*alp*ta))*e^(2*bet*ta))*e^(-2*alp*ta - 2*bet*ta)/(alp^4*bet - alp^3*bet^2 - alp^2*bet^3 + alp*bet^4) - ((alp - bet)*y_0*e^(bet*ta) + (alp*e^(alp*ta) - alp*e^(bet*ta))*th_0)^2*e^(-2*alp*ta - 2*bet*ta)/(alp - bet)^2

    ynodestates[des[index]]<-rnorm(n=1,mean=ymean,sd= sqrt(yvar))
    }

  simtrait<-ynodestates[1:n]
  return(simtrait)
  #return(c(mean(simtrait),sd(simtrait)))
  }

ABC_acc<-function(model.params,reg.params,root=root,tree=tree,summarytrait=summarytrait){
  alpha.y<-model.params[1]
  sigma.sq.y<-model.params[2]
  alpha.x<-model.params[3]
  optim.x<-model.params[4]
  sigma.sq.x<-model.params[5]

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


posterior<-run_MCMC_ABC(c(true.alpha.y,true.sigma.sq.y,true.alpha.x,true.optim.x,true.sigma.sq.x,true.b0,true.b1,true.b2),100000,root=root,tree=tree,summarytrait=summarytrait)

plot(posterior)
summary(posterior)
apply(posterior,2,FUN=mean)
plot(mean(posterior))
