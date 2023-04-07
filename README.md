# BHRM

### Introduction
BHMR approach yields conditional individual exposure-specific estimates along with an estimate of the overall effect of the mixture. While BHRM is a broad framework of statistical models, our implementation of BHRM combines :
1) a g-prior for the corresponding exposure effects to provide robust estimation of highly correlated exposures; 
2) a Bayesian selection procedure to estimate the posterior inclusion probability of individual exposures;
3) Bayesian g-computation in a potential outcome framework for estimating the overall mixture effect based on two specified exposure profiles. 

### Method
Suppose we have a model with p exposures and q covariates:
$$g(\mu_i)=\alpha+\sum_{j=1}^{p}\gamma_j\beta_jX_j+\sum_{k=1}^{q}{\delta_kU_k}+\epsilon_i,\ \ i=1,\ldots n,$$
where $X_j$ is an exposure with the corresponding estimate $\beta_j$, $\gamma_j$ is a binary variable indicating the inclusion of a specific exposure $j$ in the mixture, and $U_k$ is a covariate with corresponding effect estimates $\delta_k$. Specifically, $\boldsymbol{\gamma}=\left(\gamma_1,\ldots,\gamma_p\right)$ is a vector with element $\gamma_j\in(0,1)$  which indicates if variable $X_j$ should be included in the model $\mathcal{M}_\gamma (\gamma_j=0$ is equavalent to $\beta_j=0)$. And we chose the Beta-Binomial distribution with the function of p as the prior for $\boldsymbol{\gamma}$.

To incorporate the g-prior into the model, we interpret the regression model into a Bayesian setting under model $\mathcal{M}_\boldsymbol{\gamma}$. Assume we have a continuous outcome:

$$\boldsymbol{Y}\ |\ \alpha,\mathbf{\beta},\phi,\mathbf{\gamma}\sim N(\mathbf{1}\ \alpha+X_\mathbf{\gamma}\beta_\mathbf{\gamma},\ \phi^{-1}\boldsymbol{I}\ \ )$$

$$\mathbf{\beta}_\mathbf{\gamma}\sim\mathrm{N}\left(0,g\phi^{-1}\left({\boldsymbol{X}_\boldsymbol{\gamma}}^\prime\boldsymbol{X}_\boldsymbol{\gamma}\right)^{-1}\right)$$

where prior covariance of $\mathbf{\beta}$ is specified as the scaled version of the covariance matrix of the MLE estimator and a function of the variance of the outcome $\phi^{-1}$. We can derive the posterior distributions:

$$\mathbf{\beta}_\mathbf{\gamma}|\boldsymbol{Y},\ \alpha,\ \phi,\ \mathbf{\gamma}\sim N\ (\frac{g}{1+g}\ \widehat{\beta_\mathbf{\gamma}}, \ \phi^{-1}\frac{g}{1+g}\left({\boldsymbol{X}_\boldsymbol{\gamma}}^\prime\boldsymbol{X}_\boldsymbol{\gamma}\right)^{-1})$$

$$\alpha|\ \boldsymbol{Y},\ \mathbf{\gamma}\sim N\ (\widehat{\alpha_\mathbf{\gamma}},\ \phi^{-1}\left({\boldsymbol{X}_\boldsymbol{\gamma}}^\prime\boldsymbol{X}_\boldsymbol{\gamma}\right)^{-1})$$

where $\ \widehat{\beta_\mathbf{\gamma}}$ is the MLE estimators for $\mathcal{M}_\boldsymbol{\gamma}$. We can see that scalar g controls the conditional posterior mean shrink from the MLE estimator to prior mean zero. Also, the dispersion of the posterior covariance shrinks by the factor of $g/(1\ +\ g)$. Within this framework, the posterior inclusion probability (PIP) on the individual $\gamma_j$ is the posterior probability that the coefficient is non-zero. In this case, the Bayesian selection algorithm is combined with the shrinkage factor for g-priors to yield optimal prediction performance. For our BHMR, we adapted the semi-Bayes fixed $g$ prior that is pre-specified with a constant according to the level of the desired shrinkage.

To obtain coefficient estimations, we used the MCMC simulation method for Bayesian hierarchical models through Just another Gibbs sampler (JAGS) coding scheme. After obtaining coefficients estimations, we utilized $g$ computation to yield a single effect estimate using:(e.g., a mixtures effect) that captures the impact of one standard deviation increase in levels of all exposures simultaneously. Specifically, we use posterior predictive distributions to estimate a single mixture risk difference $(\psi_{RD})$ based on two exposure profiles, such that:
$$\psi_{RD} =  \psi_{x\ast\ =high}-\psi_{x\ast=low}.$$


### R function
We developed functions to implement the BHRM algorithm using JAGS programming framework. Users can choose the function with or without interaction effect.

To use BHRM function, input variables are defined as: 
* X: A NxP matrix of exposures for mixture analysis (on the original scale and we assume that there are no missing values)
* Y: A N-length vector for a continuous outcome
* U: A NxQ matrix of covariates (variables included in the regression model but not included in the g-estimation or model selection)
* profiles: A 2xP matrix of two counterfactual profiles of exposures for which a potential outcomes risk difference is calculated (as the exposures are assumed to be standardized, these profiles should be on the standard normal scale)
* family: a character string representing the type of outcome. Choose between "gaussian" or "binomial"
* w: weight for g-prior: w -> 0 shrink to common mean; as w -> 1 toward the maximum likelihood estimate
* n.adapt: number of iterations for the adaptation
* n.burnin: length of burn-in
* n.sample: number of random samples draws

Output variables are defined as: 
* beta: estimated coefficients for exposures effect
* beta.int: estimated coefficients for interactions effect
* gamma: PIPs for exposures
* gamma.int: PIPs for interactions
* eta.low/eta.high: contrasted profiles


For more information about how to install JAGS in R, please refer to : https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119287995.app1#:~:text=The%20following%20are%20the%20steps%20to%20install%20R%20package%20rjags,packages(%22rjags%22).



### Examples
We demonstrate the use of BHRM by following examples:

#### Model without interaction term
```{r}
library(MASS)
# run the BHRM.g R script first
# generate example data:
numInd <- 1000   # number of observations
numE <- 5   # number of dxposures
numCov <- 2   # number of covariates 
family="gaussian"   # Choose the model/outcome type from {"gaussian", "binomial"}
w <- 0.9   # Specify the g-prior weight parameter. The overall shrinkage parameter G = w/(1-w)

X <- mvrnorm(n=numInd, mu=rep(0,numE), Sigma=diag(numE))   # generate exposures data
Y <- ifelse(rep(family,numInd)=="gaussian", rnorm(numInd, 0, 1), rbinom(numInd, 1, .5))   # generate outcome data
U <- mvrnorm(n=numInd, mu=rep(0,numCov), Sigma=diag(numCov))   # generate covariates data
profiles <- rbind(rep(-0.5,numE), rep(0.5, numE))   # specify the contrasting profile (e.g., from -0.5 vs. 0.5)
n.adapt=5000; n.burnin=5000; n.sample=5000;   # specify the number of iterations for the adaptation, length of burn-in, and the number of random samples


# run BHRM analysis
results <- BHRM(X=X,Y=Y,U=U,profiles=profiles,
                family = family, w=w, selection=T,
                n.adapt=5000, n.burnin=5000, n.sample=5000)
```
<img src="https://user-images.githubusercontent.com/33040114/211767507-9d618aaa-5497-49fa-908f-c35cb8201b84.png" width="290" height="300">

#### Model with interaction term
```{r}
library(MASS)
# run the BHRM.g.interaction R script first
# generate example data:
numInd <- 1000   # number of observations
numE <- 5   # number of dxposures
numCov <- 2   # number of covariates 
family="binomial"   # Choose the model/outcome type from {"gaussian", "binomial"}
w <- 0.9    # Specify the g-prior weight parameter. The overall shrinkage parameter G = w/(1-w)

X <- mvrnorm(n=numInd, mu=rep(0,numE), Sigma=diag(numE))   # generate exposures data
Y <- ifelse(rep(family,numInd)=="gaussian", rnorm(numInd, 0, 1), rbinom(numInd, 1, .5))   # generate outcome data
U <- mvrnorm(n=numInd, mu=rep(0,numCov), Sigma=diag(numCov))   # generate covariates data
profiles <- as.data.frame(c(0,1)*matrix(1, nrow=2, ncol=numE))   # specify the contrasting profile (e.g., from -0.5 vs. 0.5)
n.adapt=5000; n.burnin=5000; n.sample=5000;   # specify the number of iterations for the adaptation, length of burn-in, and number of random samples

# run BHRM analysis with pairwise interaction terms
results <- BHRM.interaction(X=X,Y=Y,U=U,profiles=profiles,
                family = family, w=w, selection=F,
                n.adapt=5000, n.burnin=5000, n.sample=5000)
```
<img src="https://user-images.githubusercontent.com/33040114/211767549-ae6308ca-e41d-42e7-9719-d76a310e00a8.png" width="300" height="600">
