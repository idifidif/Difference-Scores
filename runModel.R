source ("functions.r")
source("CorrelationModel.r")
source("GroupsModel.r")
#Initiate the stan model - that might take some time
model<-stan_model(model_code=stanModelFull) 

#Prepare the data.
#Observed is a list of matrices, each member of which represents one group with a subjects X variables matrix. Make sure data is organized in this way. At this point the same variables should be used in all groups. The first and second variables should be O1 and O2 (such that the difference score to be calculated is O2-O1). Reliabilities - is a vector of known reliability values. Should be in the same length as the number of variables.  
dl=PrepData(Observed,Reliabilities=c(0.8,0.8,0.8)) 

#Initiate MCMC sampling
p<-sampling(object = model,data=dl,chains=3,iter=1000,warmup=200)

#Analyze sample. We suggest using either shinystan (r package), or converting the stan object into a matrix and analyzing manually. Parameters names:
#Omega is a correlation parameter - One will be usually interested in the correlation between the difference score and a third variable of interest, which will be Omega[2,3]. 2 is the indicator of the difference score, and 3 is the indicator of the first additional variable which might be added to the model (Y in the paper).
#muT - are the group means. It's a two-dimensional matrix. the first indicator is for the group, and the second for the variable. so the mean of the difference score in the second group in muT[2,2] for example.
#tau - are the group SDs. Follows the same logic as the group means
#T - are the estimates of individual participant true scores. A three-dimensional array - group X variable X subject. So for example the estimate for the true difference score of the fifth participant in the first group is T[1,2,5].
mat<-as.matrix(p)



