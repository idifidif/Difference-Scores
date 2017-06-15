stanGroupsCentral="
data{
real<lower=0,upper=1> rel1;
real<lower=0,upper=1> rel2;
int<lower=0> nSubjTot;
int<lower=0> nSubjMax;
int<lower=0> nSubj[2];
matrix[nSubj[1],2] Ogrp1;
matrix[nSubj[2],2] Ogrp2;
#int<lower=0> grp[nSubjTot];

}
parameters{
matrix [nSubj[1],2] Tgrp1;
matrix [nSubj[2],2] Tgrp2;
vector[2] muTgrp1;
vector[2] muTgrp2;
corr_matrix[2] Omegagrp1;
corr_matrix[2] Omegagrp2;
vector<lower=0>[2] taugrp1;
vector<lower=0>[2] taugrp2;
}
model{
real sigma2T2grp1;
real sigma2T2grp2;
vector[nSubj[1]] T2grp1;
vector[nSubj[2]] T2grp2;
matrix[2,2] sigma_Tgrp1;
matrix[2,2] sigma_Tgrp2;

taugrp1~cauchy(0,2.5);
taugrp2~cauchy(0,2.5);
Omegagrp1~lkj_corr(1);
Omegagrp2~lkj_corr(1);
sigma_Tgrp1=quad_form_diag(Omegagrp1,taugrp1);
sigma_Tgrp2=quad_form_diag(Omegagrp2,taugrp2);
sigma2T2grp1 = taugrp1[1]^2+taugrp1[2]^2+2*Omegagrp1[1,2]*taugrp1[1]*taugrp1[2];
sigma2T2grp2 = taugrp2[1]^2+taugrp2[2]^2+2*Omegagrp2[1,2]*taugrp2[1]*taugrp2[2];

muTgrp1~normal(0,1000);
muTgrp2~normal(0,1000);

for (s in 1:nSubj[1]){
Tgrp1[s,1:2]~multi_normal(muTgrp1,sigma_Tgrp1);
T2grp1[s]=Tgrp1[s,1]+Tgrp1[s,2];

Ogrp1[s,1]~normal(Tgrp1[s,1],sqrt((taugrp1[1]^2/rel1)-taugrp1[1]^2));
Ogrp1[s,2]~normal(T2grp1[s], sqrt((sigma2T2grp1/rel2)-sigma2T2grp1));
}    

for (s in 1:nSubj[2]){
Tgrp2[s,1:2]~multi_normal(muTgrp2,sigma_Tgrp2);
T2grp2[s]=Tgrp2[s,1]+Tgrp2[s,2];

Ogrp2[s,1]~normal(Tgrp2[s,1],sqrt((taugrp2[1]^2/rel1)-taugrp2[1]^2));
Ogrp2[s,2]~normal(T2grp2[s], sqrt((sigma2T2grp2/rel2)-sigma2T2grp2));
}

}
"


stanGroupsNonCentral="
data{
real<lower=0,upper=1> rel1;
real<lower=0,upper=1> rel2;
int<lower=0> nSubjTot;
int<lower=0> nSubjMax;
int<lower=0> nSubj[2];
matrix[nSubj[1],2] Ogrp1;
matrix[nSubj[2],2] Ogrp2;


}
parameters{
matrix [2,nSubj[1]] TSgrp1;
matrix [2,nSubj[2]] TSgrp2;
vector[2] muTgrp1;
vector[2] muTgrp2;
cholesky_factor_corr[2] omegagrp1;
cholesky_factor_corr[2] omegagrp2;
vector<lower=0>[2] taugrp1;
vector<lower=0>[2] taugrp2;
}
transformed parameters{
matrix[2,2] Omegagrp1;
matrix[2,2] Omegagrp2;
matrix [nSubj[1],2] Tgrp1;
matrix [nSubj[2],2] Tgrp2;
matrix [nSubj[1],2] muTmatgrp1;
matrix [nSubj[2],2] muTmatgrp2;
real sigma2T2grp1;
real sigma2T2grp2;
vector[nSubj[1]] T2grp1;
vector[nSubj[2]] T2grp2;

for (s in 1:nSubj[1]) for (k in 1:2) muTmatgrp1[s,k]=muTgrp1[k];
for (s in 1:nSubj[2]) for (k in 1:2) muTmatgrp2[s,k]=muTgrp2[k];
Tgrp1=muTmatgrp1+(diag_pre_multiply(taugrp1,omegagrp1)*TSgrp1)';
Tgrp2=muTmatgrp2+(diag_pre_multiply(taugrp2,omegagrp2)*TSgrp2)';
for (s in 1:nSubj[1]) T2grp1[s]=Tgrp1[s,1]+Tgrp1[s,2];
for (s in 1:nSubj[2]) T2grp2[s]=Tgrp2[s,1]+Tgrp2[s,2];

Omegagrp1=multiply_lower_tri_self_transpose(omegagrp1);
Omegagrp2=multiply_lower_tri_self_transpose(omegagrp2);
sigma2T2grp1 = taugrp1[1]^2+taugrp1[2]^2+2*Omegagrp1[1,2]*taugrp1[1]*taugrp1[2];
sigma2T2grp2 = taugrp2[1]^2+taugrp2[2]^2+2*Omegagrp2[1,2]*taugrp2[1]*taugrp2[2];
}

model{
taugrp1~cauchy(0,5);
taugrp2~cauchy(0,5);
omegagrp1~lkj_corr_cholesky(1);
omegagrp2~lkj_corr_cholesky(1);

muTgrp1~normal(0,1000);
muTgrp2~normal(0,1000);

to_vector(TSgrp1)~normal(0,1);
to_vector(TSgrp2)~normal(0,1);

for (s in 1:nSubj[1]){
Ogrp1[s,1]~normal(Tgrp1[s,1],sqrt((taugrp1[1]^2/rel1)-taugrp1[1]^2));
Ogrp1[s,2]~normal(T2grp1[s], sqrt((sigma2T2grp1/rel2)-sigma2T2grp1));
}    

for (s in 1:nSubj[2]){
Ogrp2[s,1]~normal(Tgrp2[s,1],sqrt((taugrp2[1]^2/rel1)-taugrp2[1]^2));
Ogrp2[s,2]~normal(T2grp2[s], sqrt((sigma2T2grp2/rel2)-sigma2T2grp2));
}

}
"

