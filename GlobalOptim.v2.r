
#library(BMhyb)
#GlobalOptim<-function(NegLogLike,CovRes,precision=precision,p0=p0,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2,dsm=dsm,model=model){
GlobalOptim<-function(NegLogLike,CovRes,precision=precision,p0=p0,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor= predictor,dsm=dsm,model=model){
  #names(free.parameters) <- c("alp", "bet", "sig", "xi", "b1", "b2", "alp_1", "alp_2", "alp_3", "th_1", "th_2", "sig_2")
  
  #xi is the diffusion params for the dsigma_t = """xi""" d Wt in  OUOUBM and OUBMBM model, 
  # xi is also the same to the sig_2 in dsigma_t= alp_3(th_2-sigma_t)dt + """"""sig_2"""""" sqrt(sigma_t) dWt
  #sig is the diffusion params for  dy = alpha (y-theta_t)dt + """sig""" dWt in OUBM and OUOU model
  
  #################################################
  #alp: everyone
  #bet: ouou, ououbm, ououcir
  #delta: oubmcir, ououcuir 
  #sig: oubm, ouou
  #xi: ououbm, oubmbmb, oubmcir, ououcir
  #th_1: oubmcir,ououcir
  ################################################
  
  
  #alp = alp_1,bet=alp_2,delta=alp_3
  #free.model.params<-rep(TRUE,6)  # need to check
  #names(free.model.params)<-c("alp","bet","delta","sig","xi","th_1")
  #free.pred.params<-rep(TRUE,(dim(dataset)[2]-1))
  #names(free.pred.params)<-paste("b",1:(dim(dataset)[2]-1) ,sep="")
  free.parameters<-rep(TRUE, 6)
  names(free.parameters)<-c("alp","bet","delta","sig","xi","th_1")
  
  if(model=="OUBM") {
    free.parameters[which(names(free.parameters) == "bet")]<-FALSE
    free.parameters[which(names(free.parameters) == "xi")]<-FALSE
    free.parameters[which(names(free.parameters) == "th_1")]<-FALSE
    free.parameters[which(names(free.parameters)=="delta")]<-FALSE
  }
  if(model=="OUOU") {
    free.parameters[which(names(free.parameters) == "xi")]<-FALSE
    free.parameters[which(names(free.parameters) == "th_1")]<-FALSE
    free.parameters[which(names(free.parameters)=="delta")]<-FALSE
  }
  if(model== "OUBMBM") {
    free.parameters[which(names(free.parameters) == "sig")]<-FALSE
    free.parameters[which(names(free.parameters) == "bet")]<-FALSE
    free.parameters[which(names(free.parameters) == "th_1")]<-FALSE
    free.parameters[which(names(free.parameters)=="delta")]<-FALSE
  }
  if(model== "OUOUBM") {
    free.parameters[which(names(free.parameters)=="sig")]<-FALSE
    free.parameters[which(names(free.parameters) == "th_1")]<-FALSE
    free.parameters[which(names(free.parameters)=="delta")]<-FALSE
  }
  
  if(model == "OUBMCIR") {
    free.parameters[which(names(free.parameters)=="sig")]<-FALSE
    free.parameters[which(names(free.parameters)=="bet")]<-FALSE
   
    }
  
  if(model == "OUOUCIR") {
    free.parameters[which(names(free.parameters)=="sig")]<-FALSE
  }
  
  ####CHECK START FROM HERE NEXT MEETING..............
  
  tmp_p0<-p0
  tmp_regb<-regb
 # p0<- c(2,3,4,1,1,2,5)
#  regb<-b_ini
  #best.run<-list(value=NegLogLike(p0,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm),par=p0)
  best.run<-list(value=NegLogLike(p0,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm),par=p0)
  
  print("starting likelihood")
  print(best.run$value)
  trial<-0
  impr_count<-0
  while(trial <=1){
    trial<-trial+1
    print(trial)
    #print(paste("improvement search ",trial))
    tryout<-0
    while(tryout<=3){
      tryout<-tryout+1
      #print(paste("This is the search for", trial, "th improvement where the ",tryout, "th search were trying to get the convergence estimates"))
      #print("starting points")
      #print(p0)
      if(model=="OUBM"){
        dsm<-sim.dsm.oubm(p0,predictor=predictor,x.tre=x.tre)
        new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)       
        }
      if(model=="OUOU"){
        dsm<-sim.dsm.ouou(p0,predictor=predictor,x.tre=x.tre)
        new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)       
        }
      if(model=="OUBMBM"){
        dsm<-sim.dsm.oubmbm(p0,predictor=predictor,x.tre=x.tre)
        new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)       
        }
      if(model=="OUOUBM"){
        dsm<-sim.dsm.ououbm(p0,predictor=predictor,x.tre=x.tre)
        new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)       
      }
      if(model=="OUBMCIR"){
        dsm<-sim.dsm.oubmcir(best.run$par,predictor_1=predictor_1,predictor_2=predictor_2,x.tre=x.tre)
        new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2,dsm=dsm)       
      }
      if(model=="OUOUCIR"){
        print(model)
        print(best.run$par) 
        print(sim.dsm.ououcir(best.run$par,predictor_1=predictor_1,predictor_2=predictor_2,x.tre=x.tre))
        dsm<-sim.dsm.ououcir(best.run$par,predictor_1=predictor_1,predictor_2=predictor_2,x.tre=x.tre)
        new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2,dsm=dsm)
       
        #print(new.run)
        }
      
  #    new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)       
      
      if(new.run$convergence==0){
        #print("get convergence estimates")
        #print(new.run$par)
        best.run<-new.run
        break}else{
          print("not find convergence estimates yet, start from neighbor point")
          p0<-GenerateValues(best.run$par, lower=c(0, 0, 0, 0 ,0 ,0)[which(free.parameters)], upper=c(100,100, 100, max(response)-min(response),max(response)-min(response),Inf), examined.max=10*best.run$par, examined.min=0.1*best.run$par)                 
#          c("alp","bet","delta","sig","xi","th_1")
                    #print(p0
          #print("stop here ? ")
         # c("alp", "bet", "sig", "xi", "delta, "th_1")
        }
    }#end of tryout loop
    
    #if(abs(new.run$value)+0.0001<abs(best.run$value) && new.run$value>0 ){
    #if(abs(new.run$value)<abs(best.run$value) && new.run$value>0 ){
    #  impr_count<-impr_count+1
    #  print("find an improvement")
    #  print(paste("old value is", best.run$value  ,"new improvement is" ,new.run$value))
    #  best.run<-new.run
    #  print("try to search a better improvement")
    #  p0<-GenerateValues(best.run$par, lower=c(0, 0, 0, 0, -10)[which(free.parameters)], upper=c(100,100,(range(response)[2]-range(response)[1]),(range(response)[2]-range(response)[1]),10 ), examined.max=10*best.run$par, examined.min=0.1*best.run$par)
    #  print("start point")
    #  print(p0)
    #regb<-b_est+runif(1,-1,1)
    #p0[length(best.run$par)]<-best.run$par[length(best.run$par)]
    #regb[2]<-best.run$par[length(best.run$par)]
    #print(paste("retry generate values",p0))         
    #  if(impr_count==2){print("find two improvement shall be good")
    #                    print(paste("old value is", best.run$value  ,"new improvement is" ,new.run$value))                        
    #                    ;break}
    #  }else{
    #  print("not find improvement, keep searching")
    #    print( paste("best vs new ", best.run$value,new.run$value,sep=" ") )       
    #     if(trial %%3==0){
    #       p0<-best.run$par
    #regb<-b_est+runif(1,-1,1)
    #       p0[length(best.run$par)]<-best.run$par[length(best.run$par)]
    #       #regb[2]<-best.run$par[length(best.run$par)]
    #       print("research with mls values")
    #       print(p0)
    #       }else{
    #       p0<-GenerateValues(best.run$par, lower=c(0, 0, 0, 0, -10)[which(free.parameters)], upper=c(100,100,(range(response)[2]-range(response)[1]),(range(response)[2]-range(response)[1]),10 ), examined.max=10*best.run$par, examined.min=0.1*best.run$par)   
    #       print("research with generated values")
    #       print(p0)
    #       }
    #  }
    #if(trial==2){print("search does not find any improvement")}
  }#end of trial loop 
  # minLikeIndex<-  which(Likevalue==min(Likevalue))
  
  
  
  MLEs<-best.run #so either convergence estimate or the true parameter
  print("searched likelihood value")
  print(best.run$value)
  #print(model)
  #print("the estimator values are")
  print(MLEs$par)
  
 
  
  #print(est_V)
  #we need to put differnt phyloCorrFact here
  if(model=="OUBM"){
    dsm<-sim.dsm.oubm(best.run$par,predictor=predictor,x.tre=x.tre)
    est_V<- CovRes(best.run$par,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
    }
  if(model=="OUOU"){
    dsm<-sim.dsm.ouou(best.run$par,predictor=predictor,x.tre=x.tre)
    est_V<- CovRes(best.run$par,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
  }
  if(model=="OUBMBM"){
    dsm<-sim.dsm.oubmbm(best.run$par,predictor=predictor,x.tre=x.tre)
    est_V<- CovRes(best.run$par,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
  }
  if(model=="OUOUBM"){
    dsm<-sim.dsm.ououbm(best.run$par,predictor=predictor,x.tre=x.tre)
    est_V<- CovRes(best.run$par,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
    }
  if(model=="OUBMCIR"){
    dsm<-sim.dsm.oubmcir(best.run$par,predictor_1=predictor_1,predictor_2=predictor_2,x.tre=x.tre)
    est_V<- CovRes(best.run$par,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2)
    }
  if(model=="OUOUCIR"){
    dsm<-sim.dsm.ououcir(best.run$par,predictor_1=predictor_1,predictor_2=predictor_2,x.tre=x.tre)
    est_V<- CovRes(best.run$par,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2)
    }
  est_V_inv<-pseudoinverse(est_V) 
  b_est<-pseudoinverse(t(dsm)%*%est_V_inv%*%dsm)%*%t(dsm)%*%est_V_inv%*%response
  
  result<-list(MLEs=MLEs,est_V=est_V,dsm=dsm,b_est=b_est)
  return(result)
}
