importBayes<-function(){
  setwd("/run/media/landaulab/Transcend/difscores/Bayesian")
  files<-list.files(pattern = ".Rdata")  
  Result=rep(list(list()),length(files))
  
  for (f in files){
    n<-sub(pattern = "sim",replacement = "",x = f)
    n<-as.numeric(sub(pattern = ".Rdata",replacement = "",x = n))
    load(f)
    Result[[n]]<-Res
  }
 return(Result) 
}


importSEM<-function(){
  setwd("/run/media/landaulab/Transcend/difscores/SEM")
  files<-list.files(pattern = ".Rdata")  
  ResultSEM=rep(list(list()),length(files))
  
  for (f in files){
    n<-sub(pattern = "SEM",replacement = "",x = f)
    n<-as.numeric(sub(pattern = ".Rdata",replacement = "",x = n))
    load(f)
    ResultSEM[[n]]<-ResSEM
  }
  return(ResultSEM) 
}

importOS<-function(){
  setwd("/run/media/landaulab/Transcend/difscores/OS")
  files<-list.files(pattern = ".Rdata")  
  ResultOS=rep(list(list()),length(files))
  
  for (f in files){
    n<-sub(pattern = "os",replacement = "",x = f)
    n<-as.numeric(sub(pattern = ".Rdata",replacement = "",x = n))
    load(f)
    ResultOS[[n]]<-Res
  }
  return(ResultOS) 
}


importStats<-function(){
  setwd("/run/media/landaulab/Transcend/difscores/Stats")
  files<-list.files(pattern = ".Rdata")  
  ResultStat=rep(list(list()),Nsims)
  
  for (f in files){
    n<-sub(pattern = "STAT",replacement = "",x = f)
    n<-as.numeric(sub(pattern = ".Rdata",replacement = "",x = n))
    load(f)
    ResultStat[[n]]<-ResultS
  }
  return(ResultStat) 
}