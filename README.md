# BHRM


BHMR approach yields conditional individual exposure-specific estimates along with an estimate of the overall effect of the mixture. While BHRM is a broad framework of statistical models, our implementation of BHRM combines 1) a g-prior for the corresponding exposure effects to provide robust estimation of highly correlated exposures; 2) a Bayesian selection procedure to estimate the posterior inclusion probability of individual exposures; and 3) Bayesian g-computation in a potential outcome framework for estimating the overall mixture effect based on two specified exposure profiles. Furthermore, the function incorporate a imputation procedure using normal distribution so that we can automatically impute the missings whose values are below the limit of detection (LOD).

To use BHRM function, input variables are defined as: 
* X: A NxP matrix of exposures for mixture analysis (on the original scale with NA's for individuals with below detection limit (BDL))
* Y: A N-length vector for a continuous outcome
* U: A NxQ matrix of covariates (variables included in the regression model but not included in the g-estimation)
* LOD: A P-length vector of LODs for each exposure. Individuals with missing data will have data imputed below this level of detection  
* profiles: A 2xP matrix of two counterfactual profiles of exposures for which a potential outcomes risk difference is calculated (as the expsoures are standardized wihin the funciton, these profiles should be on the standard normal scale)
* family: a character string representing the type of outcome. Choose between "gaussian" or "binomial"

An example of the function call is:

`BHRM-g(X=X.obs, Y=Y, U=U, LOD=LOD, profiles=profiles, family = "gaussian")`
