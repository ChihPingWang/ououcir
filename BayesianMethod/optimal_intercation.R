#this is with inveraction
compute.optimal.sigma.interaction<-function(opt.model.params=0, interaction=FALSE, pred.int.params = pred.int.params, pred.params.mtx=NULL ,pred.int.dataset=NULL, G = G , opt.model = opt.model){
  #given dataset, identify which preditors, interactions
  #for example r = 5 predictors, s<- 2 interactions
  #pred.int.params <- b.ini
  if(model=="OU"){
    bet<-opt.model.params[1]
  }
  prdc <- sum(stri_length(colnames(pred.int.dataset)) == 1)
  colnames(pred.int.dataset)[1:prdc]
  int.prdc <- dim(pred.int.dataset)[2] - prdc
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  inv.G<-pseudoinverse(G)
  sig1<-0
  v.pred.array<-NULL
  for(bIndex in 1:prdc){
    temp.m.pred<-c(t(one)%*%inv.G%*%pred.int.dataset[,bIndex])
    temp.m.pred<-c(temp.m.pred/(t(one)%*%inv.G%*%one))
    temp.v.pred <- t(pred.int.dataset[,bIndex] - temp.m.pred*one)%*%inv.G%*%(pred.int.dataset[,bIndex] - temp.m.pred*one)/n
    v.pred.array<-c(v.pred.array,temp.v.pred)
    for(jIndex in 1:prdc){
      temp.bij.EX.sq <- 0
      temp.bij.EX <- 0
      if(jIndex != bIndex){
        if(opt.model == "BM"){
          EX.pred <-0
          EX.sq.pred <- max(G)
          temp.bij.EX.sq <-  pred.params.mtx[bIndex,jIndex]^2*temp.v.pred*EX.sq.pred

        }
        if(opt.model == "OU"){
          x_0 <- temp.m.pred
          EX.pred <- x_0*exp(bet*max(G))
          EX.sq.pred <- (1/(2*bet))*(1 - exp(-2*bet*max(G))) + x_0^2*exp(-2*bet*max(G))# + mu*(1-exp(-bet*max(G))))^2
          temp.bij.EX.sq <-  pred.params.mtx[bIndex,jIndex]^2*temp.v.pred*EX.sq.pred
          temp.bij.EX <- pred.params.mtx[bIndex,jIndex]*EX.pred

        }#end of ou
      }
    }# jIndex
    sig1<-sig1+pred.int.params[bIndex]^2*temp.v.pred + temp.bij.EX.sq*temp.v.pred + 2*pred.int.params[bIndex]*temp.v.pred*temp.bij.EX
  }# THIS LOOP IS FOR OPTIMAL ASSUMING BM
  sig1<-sqrt(sig1)
  return(sig1)
}
