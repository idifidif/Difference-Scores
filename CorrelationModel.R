#Simple stan model centeral parametarization

T1DYCentral="
#these are the known nodes
data{ 
real<lower=0,upper=1> rel1;
real<lower=0,upper=1> rel2;
real<lower=0,upper=1> rely;
int<lower=0> nSubj;
vector[nSubj] O1;
vector[nSubj] O2;
vector[nSubj] OY;
}

#The unknown parameters:
parameters{
matrix [nSubj,3] T;
vector[3] muT;
corr_matrix[3] Omega;
vector<lower=0>[3] tau;
}

#The actual model
model{
real sigma2T2;
vector[nSubj] T2;
matrix [3,3] sigma_T;

tau~cauchy(0,5); #Hyperprior for the SDs
Omega~lkj_corr(1); #Hyperprior for the correlation matrix
muT~normal(0,1000); #Hyperprior for the mean of the three true scores distributions

sigma_T=quad_form_diag(Omega,tau); #The covariance matrix
sigma2T2 = tau[1]^2+tau[2]^2+2*Omega[1,2]*tau[1]*tau[2]; #The variance of T2.


for (s in 1:nSubj){
    T[s,1:3]~multi_normal(muT,sigma_T); #Multi-variate prior for the true scores
    T2[s]=T[s,1]+T[s,2]; #Deterministic node for T2
    
    #The likelihood distributions of the observed scores:
    O1[s]~normal(T[s,1],sqrt((tau[1]^2/rel1)-tau[1]^2));
    O2[s]~normal(T2[s], sqrt((sigma2T2/rel2)-sigma2T2));
    OY[s]~normal(T[s,3],sqrt((tau[3]^2/rely)-tau[3]^2));
  }
}
"

#(same) Stan model non-centeral parametarization. This is a more complicated definition for the same model, which is much more efficient in stan----

T1DYNonCentral="

#This part defines the known nodes of the model

data{ 
real<lower=0,upper=1> rel1; 
real<lower=0,upper=1> rel2;
real<lower=0,upper=1> rely;
int<lower=0> nSubj;
vector[nSubj] O1;
vector[nSubj] O2;
vector[nSubj] OY;
}

#This part defines the unknown parameters of the model:

parameters{
matrix [3,nSubj] TS;
vector[3] muT;
cholesky_factor_corr [3] omega;
vector<lower=0>[3] tau;
}

# More parameters (deterministic transformations of the original parameters)
transformed parameters{
matrix[3,3] Omega;
matrix [nSubj,3] T;
matrix [nSubj,3] muTmat;
vector[nSubj] T2;
real sigma2T2;
for (k in 1:3) for (s in 1:nSubj) muTmat[s,k]=muT[k];
T=muTmat+(diag_pre_multiply(tau,omega)*TS)';
for (s in 1:nSubj) T2[s]=T[s,1]+T[s,2];
Omega=multiply_lower_tri_self_transpose(omega);
sigma2T2 = tau[1]^2+tau[2]^2+2*Omega[1,2]*tau[1]*tau[2];
}

#Defining the model itself
model{
tau~cauchy(0,5); #Hyperprior for the SDs
omega~lkj_corr_cholesky(1); #Hyperprior for the correlation matrix
muT~normal(0,1000); #Hyperpriors for the means of the three variables
to_vector(TS)~normal(0,1); #Hyperpriors for the standardized true scores. This is the non-central parameterization part. Instead of drawing directly from a multivariate distribution, we draw here from standardized normal distributions and define the dependency above

#Defining the relationship between observed and true scores:
for (s in 1:nSubj){ 
O1[s]~normal(T[s,1],sqrt((tau[1]^2/rel1)-tau[1]^2));
O2[s]~normal(T2[s], sqrt((sigma2T2/rel2)-sigma2T2));
OY[s]~normal(T[s,3],sqrt((tau[3]^2/rely)-tau[3]^2));
}
}

#The covariance matrix:
generated quantities{
matrix[3,3] Sigma;
Sigma=quad_form_diag(Omega,tau);
}
"