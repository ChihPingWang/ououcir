rm(list=ls())

compute.optimal.sigma<-function(pred.params,predictor.dataset=predictor.dataset, V=V){
  r<-dim(predictor.dataset)[2]
  #free.pred.params<-rep(TRUE,r) # the maximum predcitor 
  #names(free.pred.params)<-paste("b",1:r,sep="")
  #num.pred<- dim(predictor.dataset)[2] # number of predictors
  #free.pred.params[which((names(free.pred.params) %in%  paste("b",(num.pred+1):5,sep="") ))]<- FALSE
  n<-dim(V)[1]
  one<-array(1,c(n,1))
  inv.V<-pseudoinverse(V)
  sig1<-0
#  for(bIndex in 1:length(free.pred.params)){
  for(bIndex in 1:r){
     #if(free.pred.params[bIndex]==TRUE){
      temp.m.pred<-c(t(one)%*%inv.V%*%predictor.dataset[,bIndex])
      temp.m.pred<-c(temp.m.pred/(t(one)%*%inv.V%*%one))
      temp.v.pred <- t(predictor.dataset[,bIndex] - temp.m.pred*one)%*%inv.V%*%(predictor.dataset[,bIndex] - temp.m.pred*one)/n
      sig1<-sig1+pred.params[bIndex]^2*temp.v.pred
     #}
    }
  sig1<-sqrt(sig1)
  return(sig1)
  }

size<-5
num.pred<-3
pred.params<-rexp(3, rate=5)
pred.params

predictor.dataset<-matrix(rnorm(size*num.pred),ncol=num.pred)
V<-diag(abs(rnorm(size)),c(size,size))
print(compute.optimal.sigma(pred.params, predictor.dataset = predictor.dataset, V=V))


