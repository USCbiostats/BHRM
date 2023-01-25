# Script to perform Bayesian Ridge Regression with g prior and selection for mixtures
# author: "Jingxuan He and David Conti"
# last updated: 11/29/22

library(R2jags)

# Functions
# JAGS framework that specifies the Bayesian model
BHRM.gaussian.model <- 
  "model {
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[i], prec.sigma.Y)
    mu[i] <- alpha + inprod(beta[1:P], X[i,1:P]) + inprod(delta[1:Q], U[i,1:Q])
  }
  # prior on outcome variance
  prec.sigma.Y <- 1/(sigma.Y*sigma.Y)
  sigma.Y ~ dunif(0,3)
  
  # prior on covariate effects
  for(q in 1:Q) { delta[q] ~ dnorm(0, 1.0E-06) }
  
  # prior on intercept
  alpha ~ dnorm(0, 1.0E-06)
  
  # g-prior on exposure effects
  b[1:P] ~ dmnorm(mu.beta[1:P], T[1:P, 1:P]) # multivariate normal
  for(j in 1:P) {
    mu.beta[j] <- (1-gamma[j])*prop.mu.beta[j]  # prior means
    beta[j] <- b[j]*gamma[j]                    # effect multiplied by inclusion indicator
    gamma[j] ~ dbern(pi)
    for(k in 1:P) {
      T[j,k] <- gamma[j]*gamma[k]*prec.sigma.Y*XtX[j,k]/(G) + (1-gamma[j]*gamma[k])*equals(j,k)*pow(prop.sd.beta[j],-2)
    }
  }
  
  # prior on exposure variable inclusion
  pi ~ dbeta(1,P)
  
  # semi-Bayes # w -> 0 shrink to common mean; as w -> inf toward the maximum likelihood estimate (have as an input parameter)
  G <- w/(1-w)

  # g-estimation
  eta.low <- inprod(beta[1:P], profiles[1,1:P])
  eta.high <- inprod(beta[1:P], profiles[2,1:P])
  psi <-eta.high-eta.low
  
}"


BHRM.logistic.model <- 
  "model {
  for(i in 1:N) {
    Y[i] ~ dbern(mu[i])
    logit(mu[i]) <- alpha + inprod(beta[1:P], X[i,1:P]) + inprod(delta[1:Q], U[i,1:Q])
  }

  # prior on covariate effects
  for(q in 1:Q) { delta[q] ~ dnorm(0, 1.0E-06) }
  
  # prior on intercept
  alpha ~ dnorm(0, 1.0E-06)
  
  # prior on exposure effects
  b[1:P] ~ dmnorm(mu.beta[1:P], T[1:P, 1:P])
  for(j in 1:P) {
    mu.beta[j] <- (1-gamma[j])*prop.mu.beta[j]
    beta[j] <- b[j]*gamma[j]
    gamma[j] ~ dbern(pi)
    for(k in 1:P) {
      T[j,k] <- gamma[j]*gamma[k]*XtX[j,k]/(G) + (1-gamma[j]*gamma[k])*equals(j,k)*pow(prop.sd.beta[j],-2)
    }
  }
  
  # prior on exposure variable inclusion
  pi ~ dbeta(1,P)

  # semi-Bayes # w -> 0 shrink to common mean; as w -> inf toward the maximum likelihood estimate (have as an input parameter)
  G <- w/(1-w)

  # g-estimation
  eta.low <- inprod(beta[1:P], profiles[1,1:P])
  eta.high <- inprod(beta[1:P], profiles[2,1:P])
  psi <- eta.high-eta.low
  
}"

# MCMC procedure to update the Bayesian parameters and get the estimates for the model
BHRM <- function(X=NULL, Y=NULL, U=as.matrix(rep(0, dim(X)[1])), profiles=NULL, family = "gaussian", w=0.9, n.adapt=5000, n.burnin=5000, n.sample=5000) {
  N <- length(Y)
  P <- ncol(X)
  Q <- ncol(U)
  
  univariate.results <- t(sapply(1:P, FUN=function(p) {  
    x <- as.matrix(X[,p])
    reg <- glm(Y~x+U, family=family)    # perform regression
    c.reg <- summary(reg)$coef["x",]    # select the coefficients for the exposure
  }, simplify=T))
  
  ### g prior model inputs
  prop.mu.beta <- rep(0, P)
  prop.sd.beta <- univariate.results[,"Std. Error"]   # estimate the SD using univariate results
  XtX <- t(as.matrix(X))%*%as.matrix(X) 
  
  ### run jags
  jags.model.text <- ifelse(family=="gaussian", BHRM.gaussian.model, BHRM.logistic.model)
  jdata <- list(N=N, Y=Y, X=X, U=U, P=P, Q=Q, profiles=profiles, XtX=XtX, w=w, prop.mu.beta=prop.mu.beta, prop.sd.beta=prop.sd.beta)
  var.s <- c("beta", "gamma", "eta.low", "eta.high",  "psi")
  model.fit <- jags.model(file=textConnection(jags.model.text), data=jdata, n.chains=1, n.adapt=n.adapt, quiet=T)
  update(model.fit, n.iter=n.burnin, progress.bar="none")
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=n.sample, thin=1, progress.bar="none")
  
  ### summarize results
  r <- summary(model.fit)
  BHRM.results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))
  wald = abs(BHRM.results[,"Mean"]/BHRM.results[,"SD"])
  BHRM.results$p.val = round(2*(1-pnorm(wald,0,1)), 3)
  return(BHRM.results)
  
}



