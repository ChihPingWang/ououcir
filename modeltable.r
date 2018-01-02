rm(list=ls())
library(xtable)
modelnames<-c("","OUBM","OUOU","OUBMBM","OUOUBM","OUBMCIR","OUOUCIR")
sdetype <- c("Linear","Autonomous", "Additive","Normal","Ref.")
modelproptable<-array(0,c(length(modelnames),length(sdetype)))
rownames(modelproptable)<-modelnames
colnames(modelproptable)<-sdetype
# modelproptable[,which(sdetype=="Linear")]<- "Y"
# modelproptable[,which(sdetype=="Homo")]<- "Y"
# modelproptable[,which(sdetype=="Auto")]<- "N"
# modelproptable[,which(sdetype=="Additive")]<- "Y"
# modelproptable[which(modelnames%in%c("BM","OU")),which(sdetype=="Auto")]<- "Y"
# modelproptable[,which(sdetype=="Normal")]<- "Y"
# modelproptable[which(modelnames%in%c("OUBMBM","OUOUBM","OUBMCIR","OUOUCIR")),which(sdetype=="Normal")]<- "N"
modelproptable[1,]<-"(y,t,s)"
modelproptable[2:7,which(sdetype=="Linear")]<-c(
 "(y,y,-)",
 "(y,y,-)",
 "(y,y,y)",
 "(y,y,y)",
 "(y,y,n)",
 "(y,y,n)")
 modelproptable[2:7,which(sdetype=="Autonomous")]<-c(
  "(n,y,-)",
  "(n,y,-)",
  "(n,y,y)",
  "(n,y,y)",
  "(n,y,n)",
  "(n,y,n)")
 modelproptable[2:7,which(sdetype=="Additive")]<-c(
   "(y,y,-)",
   "(y,y,-)",
   "(y,y,y)",
   "(y,y,y)",
   "(y,y,n)",
   "(y,y,n)")
modelproptable[2:7,which(sdetype=="Normal")]<-c(
     "(y,y,-)",
     "(y,y,-)",
     "(n,y,y)",
     "(n,y,y)",
     "(n,y,n)",
     "(n,y,n)")
modelproptable[2:7,which(sdetype=="Ref.")]<-c(
     "Hansen08",
     "Jhwueng14",
     "Jhwueng16",
     "Jhwueng16",
     "",
     "")
print(modelproptable)
xtable(modelproptable)
