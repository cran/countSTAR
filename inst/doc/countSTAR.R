## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----roaches------------------------------------------------------------------
data(roaches, package="countSTAR") 

# Roaches:
y = roaches$y

# Plot the PMF:
plot(0:max(y), 
     sapply(0:max(y), function(js) mean(js == y)), 
     type='h', lwd=2, main = 'PMF: Roaches Data',
     xlab = 'Roaches', ylab = 'Probability mass')

## ----roaches-covar------------------------------------------------------------
# Design matrix:
X = model.matrix( ~ roach1 + treatment + senior + log(exposure2),
                 data = roaches)

head(X)

## ----freq-lm------------------------------------------------------------------
library(countSTAR)

# EM algorithm for STAR linear regression
fit = lm_star(y ~ roach1 + treatment + senior + log(exposure2),
              data = roaches, 
              transformation = 'np')

# Fitted coefficients:
round(coef(fit), 3)

## ----conf---------------------------------------------------------------------
# Confidence interval for all coefficients
confint(fit)

## ----pval---------------------------------------------------------------------
# P-values:
round(pvals(fit), 4)

## ----predict-lm---------------------------------------------------------------
#Compute the predictive draws (just using observed points here)
y_pred = predict(fit)

## ----freq-ml------------------------------------------------------------------
#Fit STAR with random forests
suppressMessages(library(randomForest))
fit_rf = randomForest_star(y = y, X = X[,-1], # no intercept 
                           transformation = 'np')

#Fit STAR with GBM
suppressMessages(library(gbm))
fit_gbm = gbm_star(y = y, X = X[,-1], # no intercept 
                   transformation = 'np')

## ----freq-modelcomp-----------------------------------------------------------
#Look at -2*log-likelihood
-2*c(fit_rf$logLik, fit_gbm$logLik)

## ----bayes-lm, results='hide', message=FALSE, warning = FALSE-----------------
fit_blm = blm_star(y = y, X = X, 
                   transformation = 'bnp')

## ----estimates-bayes----------------------------------------------------------
# Posterior mean of each coefficient:
round(coef(fit_blm),3)

# Credible intervals for regression coefficients
ci_all_bayes = apply(fit_blm$post.beta,
      2, function(x) quantile(x, c(.025, .975)))

# Rename and print:
rownames(ci_all_bayes) = c('Lower', 'Upper')
print(t(round(ci_all_bayes, 3)))

## ----mcmcdiag, warning=FALSE, message=FALSE-----------------------------------
# MCMC diagnostics for posterior draws of the regression coefficients (excluding intercept)
plot(as.ts(fit_blm$post.beta[,-1]), 
     main = 'Trace plots', cex.lab = .75)

# (Summary of) effective sample sizes (excluding intercept)
suppressMessages(library(coda))
getEffSize(fit_blm$post.beta[,-1])

## ----modeldiag, warning=FALSE, message=FALSE----------------------------------
# Posterior predictive check using bayesplot
suppressMessages(library(bayesplot))
prop_zero = function(y) mean(y == 0)
ppc_stat(y = y, 
          yrep = fit_blm$post.pred, 
          stat = "prop_zero")

## ----bart---------------------------------------------------------------------
fit_bart = bart_star(y = y, X = X, 
                     transformation = 'np')

## ----bartppc------------------------------------------------------------------
ppc_dens_overlay(y = y, 
                 yrep = fit_bart$post.pred[1:50,])

## ----waic---------------------------------------------------------------------
waic <- c(fit_blm$WAIC, fit_bart$WAIC)
names(waic) <- c("STAR Linear Model", "BART-STAR")
print(waic)

## ----warpDLM------------------------------------------------------------------
#Visualize the data
plot(discoveries)

# Required package:
library(KFAS)

#Fit the model
warpfit = warpDLM(y = discoveries, type = "trend")

## ----warpPPC------------------------------------------------------------------
ppc_ribbon(y = as.vector(discoveries), 
           yrep = warpfit$post_pred)

