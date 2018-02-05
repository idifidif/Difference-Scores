library(lavaan)
library(MASS)
library(psych)
library(dplyr)
library(rstan)
library(coda)
library(modeest)
library(HDInterval)
library(parallel)
library(boot)
rstan_options(auto_write=TRUE)
cores<-detectCores()-1 #for parallelization



#Build structure for simulations (limiting for a positive-definite covariance matrix)-----
t=matrix(rep(NA,486*8),ncol=8)
sim=0
r1Y=0.001
for (i in c(0.5,0.9)){ #rel1
  for (j in c(0.5,0.9)){ #rel2
    for (l in c(-0.5,0.001,0.5)){ #r1D
      for (k in c(0.5,1,2)){ #varD
          for (p in c(0.001,0.4,0.8)){ #rDY
            for (n in c(50,100)){
            sim=sim+1
            sigma=matrix(c(1,l*sqrt(k),r1Y,
                           l*sqrt(k),k,p*sqrt(k),
                           r1Y,p*sqrt(k),1),
                         nrow=3,ncol=3,
                         byrow = T)
            
            rownames(sigma)<-c("T1","TD","TY")
            colnames(sigma)<-c("T1","TD","TY")
            
            #Creates a variable that sums the (signs) eigen values of the matrix, with three positive eigen values denoting a positive-denfinite matrix (which can serve as a covariance matrix. A sum of the signs of a positive-definite matrix thus should be equal to 3)
            t[sim,1]=sum(sign(eigen(sigma)$values)) 
            t[sim,2]=i
            t[sim,3]=j
            t[sim,4]=k
            t[sim,5]=l
            t[sim,6]=r1Y
            t[sim,7]=p
            t[sim,8]=n
          }}}}}}
colnames(t)<-c("Eigen","rel1","rel2","varD","r1D","r1Y","rDY","N")
t1<-data.frame(t)
t1=subset(t1,Eigen==3)  #Take only positive-definite combinations of parameters. In this case - it means all simulations


#Simulations------
#These are the actual simulations, using only combinations of parameters creating positive-definite covariance matrices. In this case this includes all possible combinations

repetitions=1000
Nsims=216
mu=c(0,0,0) #The means of T1, TD and TY are all zero
rDYsamp=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
r12samp=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
r12T=matrix(rep(NA,Nsims),ncol=1)
rDYT=matrix(rep(NA,Nsims),ncol=1)
varT2=matrix(rep(NA,Nsims),ncol=1)
relDsamp=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
corODTD=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
corResTD=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
relDformula=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
relRessamp=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
rResYsamp=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)

O1_empirical=rep(list(list()),Nsims)
O2_empirical=rep(list(list()),Nsims)
OY_empirical=rep(list(list()),Nsims)
T1Dlist=list()
rel1_empirical=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
rel2_empirical=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
rely_empirical=matrix(rep(NA,repetitions*Nsims),ncol=repetitions)
rndlist=list()


##
emp=T #Sample true scores with no sampling error (Thus sampling error exists only in the measurement error). This is conducted to isolate the effects of reliability on the results (e.g. otherwise the effect of N will incoporate both standard effects of sampling error, and its effects on the effects  of the measurement error on the final estimates and their CIs)
##

pb <- txtProgressBar(min = 0, max = Nsims, style = 3)
for (i in 1:Nsims){
  sigma=matrix(c(1.0,t1[i,5]*sqrt(t1[i,4]),t1[i,6],
                 t1[i,5]*sqrt(t1[i,4]),t1[i,4],t1[i,7]*sqrt(t1[i,4]),
                 t1[i,6],t1[i,7]*sqrt(t1[i,4]),1.0), 
               nrow=3,ncol=3,
               byrow = T) #The covariance matrix
  rownames(sigma)<-c("T1","TD","TY")
  colnames(sigma)<-c("T1","TD","TY")
  N=t1[i,8]
  rel1=t1[i,2]
  rel2=t1[i,3]
  rely=0.85
  
  
  T1D<-mvrnorm(n = N,
               mu = mu,
               Sigma = sigma,empirical = emp) #Generate random multivariate normal scores.
  
  T1Dlist[[i]]=T1D
  T1<-T1D[,1]
  TD<-T1D[,2]
  TY<-T1D[,3]
  T2<-T1+TD
  
  r12T[i]<-cor(T1,T2) #True correlation between T1 and T2. This is a simulated and not generated quantity.
  varT2[i]<-var(T2) #Variance of T2
  rDYT[i]<-cor(TD,TY) #True correlation between TD and TY - this is the quantity of interest. It should match the simulated correlation, and it will because the mvrnorm function included no sampling errors
    
  for (j in 1:repetitions){ #For each combination of parameters run many repetitions in which observed scores would be randomly generated
    O1=T1+rnorm(length(T1),0,sqrt(var(T1)/rel1-var(T1))) 
    O2=T2+rnorm(length(T2),0,sqrt(var(T2)/rel2-var(T2)))
    OY=TY+rnorm(length(TY),0,sqrt(var(TY)/rely-var(TY)))
    OD=O2-O1
    Ores=O2-O1*cor(O2,O1)*sd(O2)/sd(O1) #Residual scores
    
    rDYsamp[i,j]=cor(OD,OY)
    r12samp[i,j]=cor(O2,O1)
    rResYsamp[i,j]=cor(Ores,OY)
    relDsamp[i,j]=cor(OD,TD)^2
    corODTD[i,j]=cor(OD,TD)
    relRessamp[i,j]=cor(Ores,TD)^2
    corResTD[i,j]=cor(Ores,TD)
    
    O1_empirical[[i]][[j]]=O1
    O2_empirical[[i]][[j]]=O2
    OY_empirical[[i]][[j]]=OY
    rel1_empirical[i,j]=cor(O1,T1)^2
    rel2_empirical[i,j]=cor(O2,T2)^2
    rely_empirical[i,j]=cor(OY,TY)^2
    
  }
  
  
  setTxtProgressBar(pb, i, label=paste( round(i/Nsims*100, 0),"% done"))
}

t2<-cbind.data.frame(t1,r12T,varT2,relD=rowMeans(relDsamp),relRes=rowMeans(relRessamp),meanRDYsamp=rowMeans(rDYsamp),meanRRESYsamp=rowMeans(rResYsamp))

#Save simulation outputs to RData files (some large files...)
save(O1_empirical,O2_empirical,OY_empirical,file="Oempirical.RData")
save(T1Dlist,file="T1Dlist.RData")
save(t2,file="t2.RData")
save(rResYsamp,rDYsamp,r12samp,relDsamp,rResYsamp,corODTD,corResTD,file="correlations.RData")


#Load simulation output
load("Oempirical.Rdata")
load("t2.Rdata")
load("T1Dlist.Rdata")

#Bayesian solution----
source("CorrelationModel.R")
stanDSOOpt<-stan_model(model_code = TDTYNonCentral) #Read model
Nsims=216 #Number of scenarios
repetitions=1000 #Number of repetitions/samples per scenario

Result=rep(list(list()),Nsims)
ptm0 <- proc.time()
for (i in 1:Nsims){
  Res={}
  Warn={}
  TT=c(cor(T1Dlist[[i]][,2],T1Dlist[[i]][,3]),T1Dlist[[i]][,2]) #The true correlation
  ptm <- proc.time()
  for (rep in 1:repetitions){
    
    rel1=as.numeric(t2[i,"rel1"])
    rel2=as.numeric(t2[i,"rel2"])
    rely=0.85
    datalist=list(rel1=rel1,
                  rel2=rel2,
                  rely=rely,
                  O1=O1_empirical[[i]][[rep]],
                  O2=O2_empirical[[i]][[rep]],
                  OY=OY_empirical[[i]][[rep]],
                  nSubj=t2[i,"N"]
                  )
    w=NULL
    assign("last.warning",NULL,envir=baseenv()) #Track warnings
    ptm <- proc.time()
    stanFit<-sampling(object=stanDSOOpt,data=datalist,chains=3,iter=1000,warmup=200,thin=1)
    proc.time() - ptm
    w<-warnings() #Track warnings
    Warn[[i]]=w 
    #This algorithm checks for each model if there are errors related to poor convergence, and if so, runs another (slower) model with more strict definitions:
    if (!is.null(w)){
      assign("last.warning",NULL,envir=baseenv())
      w=NULL
      rm(stanFit)
      stanFit<-sampling(object=stanDSOOpt,data=datalist,chains=3,iter=1000,warmup=200,thin=1,
                        control=list(adapt_delta=0.99))
      w<-warnings()
      Warn[[i]]=w
      wcount=wcount+1
      print(w)
    }
    
    
    mat<-data.frame(as.matrix(stanFit)) #Convert stan object to matrix
    
    params=c("Omega.2.3.")
    for (j in 1:datalist$nSubj){params[j+1]=paste("T.",j,".2.",sep = "")} #This part extracts also the true scores of individual participnats for TD - not included in the ms
    res=data.frame(params=params) #Extract results
    res$median=sapply(seq(1,length(params),1),function(i) median(mat[,params[i]]))
    res$mean=sapply(seq(1,length(params),1),function(i) mean(mat[,params[i]]))
    res$mode=sapply(seq(1,length(params),1),function(i) asselin(mat[,params[i]]))
    res$lower=sapply(seq(1,length(params),1),function(i) hdi(mat[,params[i]])[1])
    res$upper=sapply(seq(1,length(params),1),function(i) hdi(mat[,params[i]])[2])
    res$T=TT
    
    Res[[rep]]<-res
    Result[[i]][[rep]]<-res
  }
  overallTime0=(proc.time() - ptm0)
  save(Res,file=paste("sim",i,".Rdata",sep = ""))
}
overallTime=proc.time() - ptm0
save(Result,file="AllSims.Rdata")



#SEM solution----
Nsims=216
repetitions=1000
ResultSEM=rep(list(list()),Nsims)
pb <- txtProgressBar(min = 0, max = Nsims*repetitions, style = 3)
tt=0
i=1
for (i in 1:Nsims){
  res<-list()
  for (rep in 1:repetitions){
    tt=tt+1
    O1=O1_empirical[[i]][[rep]]
    O2=O2_empirical[[i]][[rep]]
    OY=OY_empirical[[i]][[rep]]
    rel1=as.numeric(t2[i,"rel1"])
    rel2=as.numeric(t2[i,"rel2"])
    rely=0.85
    N=length(O1)
    
    SEMmodel<-paste("
                    #measurement model
                    T1=~ 1*O1
                    T2=~ 1*O2
                    TY=~ 1*OY
                    
                    O1~~",(1-rel1)*var(O1),"*O1                                                                                                                                                                                       
                    O2~~",(1-rel2)*var(O2),"*O2
                    OY~~",(1-rely)*var(OY),"*OY
                    
                    #correlations
                    T1~~TD
                    T1~~TY
                    TD~~TY
                    
                    #Difference score part
                    TD=~1*T2
                    T2~1*T1
                    T2~~0*T2
                    TD+T1~~0*T2
                    TY~~0*T2
                    
                    ",sep="",collapse="")
    
    
    df_sem=data.frame(O1,O2,OY)
    
    fit<-sem(model = SEMmodel,data=df_sem)
    fit1<-parameterestimates(fit,standardized = T) #Extract standardized parameters
    rTDTY<-as.numeric(fit1[fit1$lhs=="TY" & fit1$rhs=="TD","std.all"]) #The estimate of the correlation between TD and TY
    TDsem<-lavPredict(fit,newdata = df_sem)[,"TD"]  #This part extracts also the true scores of individual participnats for TD - not included in the ms
    
    
    
    #Bootstrapping to extract CIs on SEM parameters:
    myfun1<-function(data,indices){
      E1=(1-rel1)*var(data[indices,"O1"]) #Calculates the error variances per the bootstrapped sample
      E2=(1-rel2)*var(data[indices,"O2"])
      EY=(1-rely)*var(data[indices,"OY"])
      
      SEMmodel<-paste("
                      #measurement model
                      T1=~ 1*O1
                      T2=~ 1*O2
                      TY=~ 1*OY
                      
                      O1~~",E1,"*O1                                                                                                                                                                                          
                      O2~~",E2,"*O2
                      OY~~",EY,"*OY
                      
                      #correlations
                      T1~~TD
                      T1~~TY
                      TD~~TY
                      
                      #Difference score part
                      TD=~1*T2
                      T2~1*T1
                      T2~~0*T2
                      TD+T1~~0*T2
                      TY~~0*T2
                      
                      ",sep="",collapse="")
      df=data
      fitm<-sem(SEMmodel,data=df[indices,],
                control=list(iter.max=25),
                ,start="Mplus",check="start") #Some parameters set to speed things ups
      fit1<-parameterestimates(fitm,standardized = T) #The correlation of interest
      rTDTY<-as.numeric(fit1[fit1$lhs=="TY" & fit1$rhs=="TD","std.all"]) #This part extracts also the true scores of individual participnats for TD - not included in the ms
      Ts<-lavPredict(fitm,newdata=data)[,"TD"]
      return(cbind(rTDTY,Ts))
    }
   
    BOOT<-boot(data = df_sem,statistic = myfun1,R=100)  #Run the bootstrapping procedure (without parralelization)
    
    
    
    #This is a parallel version of the bootstrapping:
    #cl1<-makeCluster(ncores)
    #clusterEvalQ(cl1,library(lavaan))
    #clusterExport(cl1, c("rel1","rel2","rely"))
    
    #BOOT<-boot(data = df_sem,statistic = myfun1,R=500,parallel = "snow",ncpus=2,cl=cl1)
    #stopCluster(cl1)
    
    
    
    Bcor=BOOT$t[,1:N] 
    Bcor=Bcor[Bcor>=-1 & Bcor<=1] #Eliminate all bootstrapped correlations exceeding 1 or -1
    rTDTYl<-quantile(Bcor,probs = 0.025,na.rm=T) #Calculate CI
    rTDTYu<-quantile(Bcor,probs = 0.975,na.rm=T)
    
    TDboot=BOOT$t[,(N+1):(N*2)] #This part extracts also the true scores of individual participnats for TD - not included in the ms
    
    TDseml<-apply(TDboot,2,function(x) quantile(x,probs = 0.025,na.rm = T))
    TDsemu<-apply(TDboot,2,function(x) quantile(x,probs = 0.975,na.rm = T))
    
    
    res[[rep]]<-data.frame(params=c("rTDTY",paste("T",1:(N),sep="")),
                               value=c(rTDTY,TDsem),
                               lower=c(rTDTYl,TDseml),
                               upper=c(rTDTYu,TDsemu))
    save(res,file=paste("SEM",i,".Rdata",sep = ""))
    setTxtProgressBar(pb, i, label=paste( round((tt/Nsims*repetitions)*100, 0),"% done"))
 
  }
 
}



#Other solutions-----
Nsims=216
repetitions=1000

pb <- txtProgressBar(min = 0, max = Nsims, style = 3)
for (i in 1:Nsims){
  Res<-list()
  for (rep in 1:repetitions){
    O1=O1_empirical[[i]][[rep]]
    O2=O2_empirical[[i]][[rep]]
    OY=OY_empirical[[i]][[rep]]
    rel1=as.numeric(t2[i,"rel1"])
    rel2=as.numeric(t2[i,"rel2"])
    rely=0.85
    N=length(O1_empirical[[i]][[rep]])
    Res[[rep]]<-data.frame(lowerC1=rep(NA,N+1))
    Res[[rep]]$upperC1<-NA
    Res[[rep]]$observed<-NA
    Res[[rep]]$lowerC2<-NA
    Res[[rep]]$upperC2<-NA
    Res[[rep]]$Resid<-NA
    Res[[rep]]$lowerRes<-NA
    Res[[rep]]$upperRes<-NA
    
    rDY<-cor(O2-O1,OY)
    sdO1=sd(O1)
    sdO2=sd(O2)
    corO1O2=cor(O1,O2)
    relD=((sdO2/sdO1)*rel2+(sdO1/sdO2)*rel1-2*corO1O2)/ 
      ((sdO2/sdO1)+(sdO1/sdO2)-2*corO1O2) #Calculate difference scores' reliability with Equation S1
    Resid=residuals.lm(lm(O2~O1)) #Residual scores
    rResY=cor(Resid,OY) #correlation between residual scores and observed scores
    
    
    relD=max(relD,0.00001) #If reliability of difference scores is lower than 0 - conver to 0.00001
    Res[[rep]]$relD=rep(relD,N+1)
    
    Res[[rep]][1,"observed"]<-rDY
    Res[[rep]][1,"lowerC1"]<-r.con(rDY,n = N)[1]      
    Res[[rep]][1,"upperC1"]<-r.con(rDY,n = N)[2]
    Res[[rep]][1,"C2"]<-max(min(rDY/sqrt(relD*rely),1),-1) #Spearman correction
    Res[[rep]][1,"lowerC2"]<-max(min(Res[[rep]][1,"lowerC1"]/sqrt(relD*rely),1),-1)  #Spearman correction CI
    Res[[rep]][1,"upperC2"]<-max(min(Res[[rep]][1,"upperC1"]/sqrt(relD*rely),1),-1)
    Res[[rep]][1,"Resid"]<-rResY #Correlation based on residuals solution
    Res[[rep]][1,"lowerRes"]<-r.con(rResY,n = N)[1] #CI for residuals solution
    Res[[rep]][1,"upperRes"]<-r.con(rResY,n = N)[2]
   }
 
    
  
  save(Res,file=paste("os",i,".Rdata",sep = ""))
  setTxtProgressBar(pb, i, label=paste( round(i/Nsims*100, 0),"% done"))
}




#Load solutions' results
source("importing.R")
Result<-importBayes() #These functions are included in a separate file
ResultSEM<-importSEM()
ResultOS<-importOS()
load("t2.Rdata")


#Calculate stats-----
Nsims=216
repetitions=1000
ResultStat=rep(list(list()),Nsims)


pb <- txtProgressBar(min = 0, max = Nsims, style = 3)
for (i in 1:Nsims){
  ResultS<-list()
  for (rep in 1:repetitions){
    N=nrow(Result[[i]][[rep]])-1
  #  ResultSEM[[i]][[rep]][1,"value"]<-max(0,min(ResultSEM[[i]][[rep]][1,"value"],1))
    
    if(ResultSEM[[i]][[rep]]$value[1]>1) {ResultSEM[[i]][[rep]]$value[1]=
      (ResultSEM[[i]][[rep]]$upper[1]+ResultSEM[[i]][[rep]]$lower[1])/2}
    if(ResultSEM[[i]][[rep]]$value[1]<(-1)) {ResultSEM[[i]][[rep]]$value[1]=
      (ResultSEM[[i]][[rep]]$upper[1]+ResultSEM[[i]][[rep]]$lower[1])/2}
    #if SEM observed correlation is above 1 or below -1, use the median of the CI instead
    
    #Bayesian bias
    ResultS[[rep]]<-data.frame(ModeBias=rep(NA,N+1))
    ResultS[[rep]]$ModeBias=Result[[i]][[rep]]$mode-Result[[i]][[rep]]$T
    ResultS[[rep]]$MedianBias=Result[[i]][[rep]]$median-Result[[i]][[rep]]$T
    ResultS[[rep]]$MeanBias=Result[[i]][[rep]]$mean-Result[[i]][[rep]]$T
    ResultS[[rep]]$SEMBias=ResultSEM[[i]][[rep]]$value-Result[[i]][[rep]]$T #SEM bias
    ResultS[[rep]]$C2Bias=ResultOS[[i]][[rep]]$C2-Result[[i]][[rep]]$T #Spearman correcction bias
    ResultS[[rep]]$observedBias=ResultOS[[i]][[rep]]$observed-Result[[i]][[rep]]$T #Observed bias
    ResultS[[rep]]$ResidBias=ResultOS[[i]][[rep]]$Resid-Result[[i]][[rep]]$T #Residual solution bias
    
    #Same for RMSE (square-rooting is done after averaing in the next section)
    ResultS[[rep]]$ModeRMSE=(Result[[i]][[rep]]$mode-Result[[i]][[rep]]$T)^2
    ResultS[[rep]]$MedianRMSE=(Result[[i]][[rep]]$median-Result[[i]][[rep]]$T)^2
    ResultS[[rep]]$MeanRMSE=(Result[[i]][[rep]]$mean-Result[[i]][[rep]]$T)^2
    ResultS[[rep]]$SEMRMSE=(ResultSEM[[i]][[rep]]$value-Result[[i]][[rep]]$T)^2
    ResultS[[rep]]$C2RMSE=(ResultOS[[i]][[rep]]$C2-Result[[i]][[rep]]$T)^2
    ResultS[[rep]]$observedRMSE=(ResultOS[[i]][[rep]]$observed-Result[[i]][[rep]]$T)^2
    ResultS[[rep]]$ResidRMSE=(ResultOS[[i]][[rep]]$Resid-Result[[i]][[rep]]$T)^2
    
    
    #Coverage calculation - whether CI includes the true correlation - Bayesian
    ResultS[[rep]]$BCoverage=(Result[[i]][[rep]]$T<=Result[[i]][[rep]]$upper &
                                    Result[[i]][[rep]]$T>=Result[[i]][[rep]]$lower)
    
    #Coverage calculation - observed
    ResultS[[rep]]$C1Coverage=(Result[[i]][[rep]]$T<=ResultOS[[i]][[rep]]$upperC1 &
                                     Result[[i]][[rep]]$T>=ResultOS[[i]][[rep]]$lowerC1)
    
    #Coverage calculation Spearman correction
    ResultS[[rep]]$C2Coverage=(Result[[i]][[rep]]$T<=ResultOS[[i]][[rep]]$upperC2 &
                                    Result[[i]][[rep]]$T>=ResultOS[[i]][[rep]]$lowerC2)
    
    #Coverage calculation SEM
    ResultS[[rep]]$SEMCoverage=(Result[[i]][[rep]]$T<=ResultSEM[[i]][[rep]]$upper &
                                         Result[[i]][[rep]]$T>=ResultSEM[[i]][[rep]]$lower)
    
    #Coverage calculation Residual solution
    ResultS[[rep]]$residCoverage=(Result[[i]][[rep]]$T<=ResultOS[[i]][[rep]]$upperRes &
                                  Result[[i]][[rep]]$T>=ResultOS[[i]][[rep]]$lowerRes)
    
    
    #CI width - Bayesian:
    ResultS[[rep]]$BWidth=(Result[[i]][[rep]]$upper -Result[[i]][[rep]]$lower) 
    #CI width - Obseved:
    ResultS[[rep]]$C1Width=(ResultOS[[i]][[rep]]$upperC1 -ResultOS[[i]][[rep]]$lowerC1)
    #CI width - Spearman:
    ResultS[[rep]]$C2Width=(ResultOS[[i]][[rep]]$upperC2 -ResultOS[[i]][[rep]]$lowerC2)
    #CI width - SEM:
    ResultS[[rep]]$SEMWidth=(ResultSEM[[i]][[rep]]$upper -ResultSEM[[i]][[rep]]$lower)
    #CI width - Residuals:
    ResultS[[rep]]$ResidWidth=(ResultOS[[i]][[rep]]$upperRes -ResultOS[[i]][[rep]]$lowerRes)
    

  }
  save(ResultS,file=paste("STAT",i,".Rdata",sep = ""))
  setTxtProgressBar(pb, i, label=paste( round(i/Nsims*100, 0),"% done"))
  
}


#Summarize-----

Nsims=216
repetitions=1000
ResultStat<-importStats() #Import results into a single file
load("t2.Rdata")

#Average results across repetitions/samples for each scenario
Cor<-t2[1:Nsims,2:ncol(t2)]
rep=1
for (i in 1:Nsims){
  for (j in 1:ncol(ResultStat1[[i]][[rep]])){
    #   N=nrow(Result[[i]][[rep]])-1
    name=colnames(ResultStat1[[i]][[rep]])[j]
    value=mean(unlist(lapply(ResultStat1[[i]],function(x) mean(x[1,name]))))
    Cor[[name]][i]<-value
  }
}



Cor2<-Cor


#Take square root for RMSE calculations
Cor[,grepl("RMSE",colnames(Cor),fixed=T)]<-sqrt(Cor[,grepl("RMSE",colnames(Cor),fixed=T)]) 
Scores[,grepl("RMSE",colnames(Cor),fixed=T)]<-sqrt(Scores[,grepl("RMSE",colnames(Cor),fixed=T)])



#Calculate average results across different scenarios
library(dplyr)
Gvars=c("rel1","rel2","varD","r1D","rDY","N")
Cor1<-Cor
for (v in Gvars){
tt<-Cor1%>%
  group_by_(.dots=Gvars[Gvars!=v])%>%
  summarise_all(funs(mean))
  tt[v]<-"Average"
  Cor1<-rbind.data.frame(Cor1,tt)
}

write.csv(Cor1,"corSEMsec.csv")


#Plots for 10th, 50th, and 90th percentiles
corODOY=matrix(rep(NA,216*1000),nrow=216)
for (sim in 1:216){
  for (rep in 1:1000){
  corODOY[sim,rep]=cor(O2_empirical[[sim]][[rep]]-O1_empirical[[sim]][[rep]],OY_empirical[[sim]][[rep]])
  }
}


plots<-data.frame(sim=rep(1:216,each=3),quant=rep(c("10","50","90"),216),rTDTY=NA,
                  rODOY=NA,rODOYl=NA,rODOYu=NA,
                  Bayes=NA,Bayesl=NA,Bayesu=NA,
                  SEM=NA,SEMl=NA,SEMu=NA,
                  Spearman=NA,Spearmanl=NA,Spearmanu=NA,
                  Resid=NA,Residl=NA,Residl=NA)

at50={}
at10={}
at90={}

for (sim in 1:216){
  at50[sim]=which.min(abs(corODOY[sim,]-quantile(x = corODOY[sim,],probs = 0.5)))
  at10[sim]=which.min(abs(corODOY[sim,]-quantile(x = corODOY[sim,],probs = 0.1)))
  at90[sim]=which.min(abs(corODOY[sim,]-quantile(x = corODOY[sim,],probs = 0.9)))
  N=length(O2_empirical[[sim]][[1]])
  plots[plots$sim==sim & plots$quant=="10","rODOY"]<-cor(O2_empirical[[sim]][[at10[sim]]]-O1_empirical[[sim]][[at10[sim]]],
                                                         OY_empirical[[sim]][[at10[sim]]])
  plots[plots$sim==sim & plots$quant=="10","rODOYl"]<-r.con(plots[plots$sim==sim & plots$quant=="10","rODOY"],n=N)[1]
  plots[plots$sim==sim & plots$quant=="10","rODOYu"]<-r.con(plots[plots$sim==sim & plots$quant=="10","rODOY"],n=N)[2]
  
  plots[plots$sim==sim & plots$quant=="50","rODOY"]<-cor(O2_empirical[[sim]][[at50[sim]]]-O1_empirical[[sim]][[at50[sim]]],
                                                         OY_empirical[[sim]][[at50[sim]]])
  plots[plots$sim==sim & plots$quant=="50","rODOYl"]<-r.con(plots[plots$sim==sim & plots$quant=="50","rODOY"],n=N)[1]
  plots[plots$sim==sim & plots$quant=="50","rODOYu"]<-r.con(plots[plots$sim==sim & plots$quant=="50","rODOY"],n=N)[2]
  
  plots[plots$sim==sim & plots$quant=="90","rODOY"]<-cor(O2_empirical[[sim]][[at90[sim]]]-O1_empirical[[sim]][[at90[sim]]],
                                                         OY_empirical[[sim]][[at90[sim]]])
  plots[plots$sim==sim & plots$quant=="90","rODOYl"]<-r.con(plots[plots$sim==sim & plots$quant=="90","rODOY"],n=N)[1]
  plots[plots$sim==sim & plots$quant=="90","rODOYu"]<-r.con(plots[plots$sim==sim & plots$quant=="90","rODOY"],n=N)[2]
}

#Result<-importBayes()

for (sim in 1:216){
  N=length(O2_empirical[[sim]][[1]])
  plots[plots$sim==sim & plots$quant=="10","rTDTY"]<-Result[[sim]][[at10[sim]]]$T[1]
  plots[plots$sim==sim & plots$quant=="10","Bayes"]<-Result[[sim]][[at10[sim]]]$median[1]
  plots[plots$sim==sim & plots$quant=="10","Bayesl"]<-Result[[sim]][[at10[sim]]]$lower[1]
  plots[plots$sim==sim & plots$quant=="10","Bayesu"]<-Result[[sim]][[at10[sim]]]$upper[1]
  
  plots[plots$sim==sim & plots$quant=="50","rTDTY"]<-Result[[sim]][[at50[sim]]]$T[1]
  plots[plots$sim==sim & plots$quant=="50","Bayes"]<-Result[[sim]][[at50[sim]]]$median[1]
  plots[plots$sim==sim & plots$quant=="50","Bayesl"]<-Result[[sim]][[at50[sim]]]$lower[1]
  plots[plots$sim==sim & plots$quant=="50","Bayesu"]<-Result[[sim]][[at50[sim]]]$upper[1]
  
  plots[plots$sim==sim & plots$quant=="90","rTDTY"]<-Result[[sim]][[at90[sim]]]$T[1]
  plots[plots$sim==sim & plots$quant=="90","Bayes"]<-Result[[sim]][[at90[sim]]]$median[1]
  plots[plots$sim==sim & plots$quant=="90","Bayesl"]<-Result[[sim]][[at90[sim]]]$lower[1]
  plots[plots$sim==sim & plots$quant=="90","Bayesu"]<-Result[[sim]][[at90[sim]]]$upper[1]
}

#ResultSEM<-importSEM()

#ResultSEM=ResultSEM2
#sim=1



for (sim in 1:216){
  N=length(O2_empirical[[sim]][[1]])
 
  plots[plots$sim==sim & plots$quant=="10","SEM"]<-ResultSEM[[sim]][[at10[sim]]]$value[1]
  plots[plots$sim==sim & plots$quant=="10","SEMl"]<-ResultSEM[[sim]][[at10[sim]]]$lower[1]
  plots[plots$sim==sim & plots$quant=="10","SEMu"]<-ResultSEM[[sim]][[at10[sim]]]$upper[1]
  
  
  plots[plots$sim==sim & plots$quant=="50","SEM"]<-ResultSEM[[sim]][[at50[sim]]]$value[1]
  plots[plots$sim==sim & plots$quant=="50","SEMl"]<-ResultSEM[[sim]][[at50[sim]]]$lower[1]
  plots[plots$sim==sim & plots$quant=="50","SEMu"]<-ResultSEM[[sim]][[at50[sim]]]$upper[1]
  
  
  plots[plots$sim==sim & plots$quant=="90","SEM"]<-ResultSEM[[sim]][[at90[sim]]]$value[1]
  plots[plots$sim==sim & plots$quant=="90","SEMl"]<-ResultSEM[[sim]][[at90[sim]]]$lower[1]
  plots[plots$sim==sim & plots$quant=="90","SEMu"]<-ResultSEM[[sim]][[at90[sim]]]$upper[1]
}


#ResultOS<-importOS()
for (sim in 1:216){
  N=length(O2_empirical[[sim]][[1]])
 
  plots[plots$sim==sim & plots$quant=="10","Spearman"]<-ResultOS[[sim]][[at10[sim]]]$C2[1]
  plots[plots$sim==sim & plots$quant=="10","Spearmanl"]<-ResultOS[[sim]][[at10[sim]]]$lowerC2[1]
  plots[plots$sim==sim & plots$quant=="10","Spearmanu"]<-ResultOS[[sim]][[at10[sim]]]$upperC2[1]
  
  
  plots[plots$sim==sim & plots$quant=="50","Spearman"]<-ResultOS[[sim]][[at50[sim]]]$C2[1]
  plots[plots$sim==sim & plots$quant=="50","Spearmanl"]<-ResultOS[[sim]][[at50[sim]]]$lowerC2[1]
  plots[plots$sim==sim & plots$quant=="50","Spearmanu"]<-ResultOS[[sim]][[at50[sim]]]$upperC2[1]
  
  
  plots[plots$sim==sim & plots$quant=="90","Spearman"]<-ResultOS[[sim]][[at90[sim]]]$C2[1]
  plots[plots$sim==sim & plots$quant=="90","Spearmanl"]<-ResultOS[[sim]][[at90[sim]]]$lowerC2[1]
  plots[plots$sim==sim & plots$quant=="90","Spearmanu"]<-ResultOS[[sim]][[at90[sim]]]$upperC2[1]
  
  
  plots[plots$sim==sim & plots$quant=="10","Resid"]<-ResultOS[[sim]][[at10[sim]]]$Resid[1]
  plots[plots$sim==sim & plots$quant=="10","Residl"]<-ResultOS[[sim]][[at10[sim]]]$lowerRes[1]
  plots[plots$sim==sim & plots$quant=="10","Residu"]<-ResultOS[[sim]][[at10[sim]]]$upperRes[1]
  
  
  plots[plots$sim==sim & plots$quant=="50","Resid"]<-ResultOS[[sim]][[at50[sim]]]$Resid[1]
  plots[plots$sim==sim & plots$quant=="50","Residl"]<-ResultOS[[sim]][[at50[sim]]]$lowerRes[1]
  plots[plots$sim==sim & plots$quant=="50","Residu"]<-ResultOS[[sim]][[at50[sim]]]$upperRes[1]
  
  
  plots[plots$sim==sim & plots$quant=="90","Resid"]<-ResultOS[[sim]][[at90[sim]]]$Resid[1]
  plots[plots$sim==sim & plots$quant=="90","Residl"]<-ResultOS[[sim]][[at90[sim]]]$lowerRes[1]
  plots[plots$sim==sim & plots$quant=="90","Residu"]<-ResultOS[[sim]][[at90[sim]]]$upperRes[1]
}

write.csv(plots,"plots.csv")

