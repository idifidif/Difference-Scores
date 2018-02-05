# Difference-Scores
Code for Bayesian modelling of difference scores. 
See Fradkin, I., Siegelman, N & Huppert, J. D. (Submitted). Title


# File list for the paper
-fullCode.r - includes the full code for the simulations, and analyses.

- CorrelationModel.r - includes only the STAN models focused on the analysis of the correlation between difference scores and a third variable of interest (i.e. TY). Central parameterization is given to enhance the understanding of the model structure. Non-central parameterization is more efficient. This is used by fullCode.r

- importing.r - several additional functions used in the fullcode.r for importing saved results etc.

# Files for flexible models and analysis
This includes additional files for Bayesian modeling of difference scores, which could be used on continous or categorical TY, and on an unlimited number of such variables.

- runModel.r - the main file to run the model specified in 'functions.r'. Includes instructions for data analysis.

- functions.r - includes code for a general STAN model which can be used to analysis of both correlation with additional variables, and group differences. Also includes a function to convert a datalist of the data one wants to analyze (each member is a group, and contains a matrix of subject X variable) into an array which can be analyzed by the STAN model supplied.

# Additional files
- GroupsModel.r - includes STAN models focused on the analysis of the relationship between difference scores and a dichotomous third variable of interest - group. That is: group differences in difference scores. Central parameterization is given to enhance the understanding of the model structure. Non-central parameterization is more efficient.



