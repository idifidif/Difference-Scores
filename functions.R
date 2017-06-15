library(rstan)

#STAN------
stanModelFull="
  data{
    int J; #Number of variables
    int K; #Number of groups
    vector<lower=0,upper=1>[J] rel; #Known reliability values of variables
    int<lower=0> nSubj[K+1]; #Vector of numbers of subjects per group
    int<lower=0> nSubjMax;
    matrix[nSubjMax,J] O[K]; #Creates a GroupXSubjectXVariable array
    
  }

parameters{
  matrix [J,nSubjMax] TS[K];
  vector[3] muT[K]; # Group means (T1,TD, and all additional variables)
  cholesky_factor_corr [3] omega[K];
  vector<lower=0>[3] tau[K]; #Group SDs (T1,TD, and all additional variables)
}

transformed parameters{
  matrix[3,3] Omega[K]; #Array of correlation matrices. One per group. 
  matrix [nSubjMax,3] T[K]; #Array of inidivual participants' estimated true scores - GroupXSubjectXVariable
  matrix [nSubjMax,3] muTmat[K];
  vector[nSubjMax] T2[K];
  real sigma2T2[K];
  
for (k in 1:K) {
  for (s in 1:nSubj[k])  for (i in 1:3) muTmat[k,s,i]=muT[k,i];
  T[k,,]=muTmat[k,,]+(diag_pre_multiply(tau[k,],omega[k,,])*TS[k,,])';
  for (s in 1:nSubj[k]) T2[k,s]=T[k,s,1]+T[k,s,2];
  Omega[k,,]=multiply_lower_tri_self_transpose(omega[k,,]);
  sigma2T2[k] = tau[k,1]^2+tau[k,2]^2+2*Omega[k,1,2]*tau[k,1]*tau[k,2];
}
}

model{
  for (k in 1:K){
    tau[k,]~cauchy(0,5);
    omega[k,,]~lkj_corr_cholesky(1);
    muT[k,]~normal(0,1000);
    to_vector(TS[k,,])~normal(0,1);
    for (s in 1:nSubj[k]){
    O[k,s,1]~normal(T[k,s,1],sqrt((tau[k,1]^2/rel[1])-tau[k,1]^2));
    O[k,s,2]~normal(T2[k,s], sqrt((sigma2T2[k]/rel[2])-sigma2T2[k]));
    if (J>2) for (j in 3:J) O[k,s,j]~normal(T[k,s,j],sqrt((tau[k,j]^2/rel[j])-tau[k,j]^2));
    }
  }
}

generated quantities{
if (K>2){
matrix [K,J-2] rDY;
for (k in 1:K) for (j in 3:J) rDY[k,j]=Omega[k,2,j];
}
}

"



PrepData<-function(Observed,Reliabilities){
  K=length(Observed)
  J=ncol(Observed[[1]])
  nSubj={}
  for (k in 1:K) nSubj[k]=nrow(Observed[[k]])
  nSubjMax=max(nSubj)
  Otemp={}
  for (k in 1:K) Otemp=c(Otemp,as.vector(Observed[[k]]))
  O=array(Otemp,dim=c(nSubjMax,J,K))
  OO<-aperm(O,c(3,1,2))
  datalist=list(K=K,
                J=J,
                rel=pmin(Reliabilities,0.999),
                nSubj=c(nSubj,0),
                nSubjMax=nSubjMax,
                O=OO
  )
  return (datalist)
  
}

