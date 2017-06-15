# Difference-Scores
Code for Bayesian modelling of difference scores. 
See Fradkin, I., Siegelman, N & Huppert, J. D. (Submitted). Novel solutions for (very) old problems: using Bayesian modeling for the analysis of difference scores.


We update this repository constantly, to add new models and update old ones. So keep updated. For additional information - idifidif@gmail.com

# File list:
- CorrelationModel.r - includes STAN models focused on the analysis of the correlation between difference scores and a third variable of interest (i.e. TY). Central parameterization is given to enhance the understanding of the model structure. Non-central parameterization is more efficient.

- GroupsModel.r - includes STAN models focused on the analysis of the relationship between difference scores and a dichotomous third variable of interest - group. That is: group differences in difference scores. Central parameterization is given to enhance the understanding of the model structure. Non-central parameterization is more efficient.

- functions.r - includes code for a general STAN model which can be used to analysis of both correlation with additional variables, and group differences. Also includes a function to convert a datalist of the data one wants to analyze (each member is a group, and contains a matrix of subject X variable) into an array which can be analyzed by the STAN model supplied.


- runModel.r - the main file to run the model specified in 'functions.r'. Includes instructions for data analysis.
