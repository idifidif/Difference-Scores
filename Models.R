#STAN------

stanString="
  data{
    real<lower=0,upper=1> rel1;
    real<lower=0,upper=1> rel2;
    real<lower=0,upper=1> rely;
    int<lower=0> nSubj;
    vector[nSubj] O1;
    vector[nSubj] O2;
    vector[nSubj] OY;
  }
parameters{
  matrix [nSubj,3] T;
  vector[3] muT;
  corr_matrix[3] Omega;
  vector<lower=0>[3] tau;
}
model{
  real sigma2T2;
  vector[nSubj] T2;
  matrix [3,3] sigma_T;
  
  tau~cauchy(0,2.5);
  Omega~lkj_corr(1);
  
  muT[1]~normal(0,1000);
  muT[2]~normal(0,1000);
  muT[3]~normal(0,1000);

  sigma_T=quad_form_diag(Omega,tau);
  sigma2T2 = tau[1]^2+tau[2]^2+2*Omega[1,2]*tau[1]*tau[2];

  for (s in 1:nSubj){
    T[s,1:3]~multi_normal(muT,sigma_T);
    T2[s]=T[s,1]+T[s,2];
    O1[s]~normal(T[s,1],sqrt((tau[1]^2/rel1)-tau[1]^2));
    O2[s]~normal(T2[s], sqrt((sigma2T2/rel2)-sigma2T2));
    OY[s]~normal(T[s,3],sqrt((tau[3]^2/rely)-tau[3]^2));
  }
 

}
"


stanModelFull="
  data{
    int J;
    int K;
    vector<lower=0,upper=1>[J] rel;
    int<lower=0> nSubj[K+1];
    int<lower=0> nSubjMax;
    matrix[nSubjMax,J] O[K];
    
  }

parameters{
  matrix [J,nSubjMax] TS[K];
  vector[3] muT[K];
  cholesky_factor_corr [3] omega[K];
  vector<lower=0>[3] tau[K];
}

transformed parameters{
  matrix[3,3] Omega[K];
  matrix [nSubjMax,3] T[K];
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
generated quantities{
matrix [K,J-2] rDY;
for (k in 1:K) for (j in 3:J) rDY[k,j]=Omega[k,2,j]
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
"

library(rstan)
library(MASS)
model<-stan_model(model_code=stanStringOpt)


sigma=matrix(c(1,0.3,0.2,
               0.3,1,0.7,
               0.2,0.7,1),nrow=3,ncol=3)
mu=c(0,0,0)

sigma
T=mvrnorm(50,mu,sigma,empirical=T)
T2=T[,1]+T[,2]
rel=c(0.8,0.8,0.8)
O1<-rnorm(50,T[,1],sqrt(var(T[,1])/rel[1]-var(T[,1])))
O2<-rnorm(50,T2,sqrt(var(T2)/rel[2]-var(T2)))
OY<-rnorm(50,T[,3],sqrt(var(T[,3])/rel[3]-var(T[,3])))

mu2=c(10,20,30)
Tgrp2=mvrnorm(50,mu2,sigma,empirical=T)
T2grp2=Tgrp2[,1]+Tgrp2[,2]
O1grp2<-rnorm(50,Tgrp2[,1],sqrt(var(Tgrp2[,1])/rel[1]-var(Tgrp2[,1])))
O2grp2<-rnorm(50,T2grp2,sqrt(var(T2grp2)/rel[2]-var(T2grp2)))
OYgrp2<-rnorm(50,Tgrp2[,3],sqrt(var(Tgrp2[,3])/rel[3]-var(Tgrp2[,3])))



O=array(c(O1,O2,OY,O1grp2,O2grp2,OYgrp2),dim=c(50,3,2))
OO<-aperm(O,c(3,1,2))


datalist=list(K=2,
              J=3,
              rel=c(0.99,0.99,0.99),
              nSubj=c(50,50),
              nSubjMax=max(c(50,50)),
              O=OO
              )

p<-sampling(object = model,data=datalist,chains=3,iter=1000,warmup=200)
mat<-as.matrix(p)

colnames(mat)

plot(density(mat[,"Omega[1,2,3]"]))
cor(O2-O1,OY)

plot(density(mat[,"Omega[2,2,3]"]))
cor(O2grp2-O1grp2,OYgrp2)

plot(density(mat[,"muT[2,1]"]))
