## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----roaches------------------------------------------------------------------
data(roaches, package="countSTAR") 

# Roaches:
y = roaches$y

# Function to plot the point mass function:
stickplot = function(y, ...){
  js = 0:max(y); 
  plot(js, 
       sapply(js, function(js) mean(js == y)), 
       type='h', lwd=2, ...)
}
stickplot(y, main = 'PMF: Roaches Data',
          xlab = 'Roaches', ylab = 'Probability mass')

## ----freq-lm------------------------------------------------------------------
library(countSTAR)

# Select a transformation:
transformation = 'np' # Estimated transformation using empirical CDF

# EM algorithm for STAR (using the log-link)
fit_em = lm_star(y ~ roach1 + treatment + senior + log(exposure2),
                 data = roaches, transformation = transformation)


# Dimensions:
n = nrow(fit_em$X); p = ncol(fit_em$X)

# Fitted coefficients:
round(coef(fit_em), 3)

## ----conf---------------------------------------------------------------------
# Confidence interval for all coefficients
confint(fit_em)

## ----pval---------------------------------------------------------------------
# P-values:
print(pvals(fit_em))

## ----predict-lm---------------------------------------------------------------
#Compute the predictive draws (just using observed points here)
y_pred = predict(fit_em)

## ----freq-ml------------------------------------------------------------------
# Select a transformation:
transformation = 'np' # Estimated transformation using empirical CDF

# Construct data matrix
y = roaches$y
X = roaches[, c("roach1", "treatment", "senior", "exposure2")]

#Fit STAR with random forests
fit_rf = randomForest_star(y, X, transformation = transformation)

#Fit STAR with GBM
fit_gbm = gbm_star(y, X, transformation = transformation)

## ----freq-modelcomp-----------------------------------------------------------
#Look at -2*log-likelihood
print(-2*c(fit_rf$logLik, fit_gbm$logLik))

## ----bayes-lm, results='hide', message=FALSE----------------------------------
X = model.matrix(y ~ roach1 + treatment + senior + log(exposure2),
                 data = roaches)

# Dimensions:
n = nrow(X); p = ncol(X)

fit_blm = blm_star(y = y, X=X, transformation = 'np')

## ----estimates-bayes----------------------------------------------------------
# Posterior mean of each coefficient:
round(coef(fit_blm),3)

# Credible intervals for regression coefficients
ci_all_bayes = apply(fit_blm$post.beta,
      2, function(x) quantile(x, c(.025, .975)))

# Rename and print:
rownames(ci_all_bayes) = c('Lower', 'Upper')
print(t(round(ci_all_bayes, 3)))

## ----diag, warning=FALSE, message=FALSE---------------------------------------
# MCMC diagnostics for posterior draws of the regression coefficients
plot(as.ts(fit_blm$post.beta), main = 'Trace plots', cex.lab = .75)

# (Summary of) effective sample sizes across coefficients:
getEffSize(fit_blm$post.beta)

# Posterior predictive check using bayesplot
suppressMessages(library(bayesplot))
prop_zero <- function(y) mean(y == 0)
(ppc_stat(y=roaches$y, yrep=fit_blm$post.pred, stat = "prop_zero"))

## ----bart---------------------------------------------------------------------
#Get the model matrix of predictors (no intercept necessary)
X = model.matrix(y ~ -1 + roach1 + treatment + senior + exposure2,
                 data = roaches)

fit_bart = bart_star(y = y, X=X, transformation = 'np')

## ----bartppc------------------------------------------------------------------
ppc_dens_overlay(y=roaches$y, yrep=fit_bart$post.pred[1:50,])

## ----waic---------------------------------------------------------------------
waic <- c(fit_blm$WAIC, fit_bart$WAIC)
names(waic) <- c("STAR w/ Linear Model", "STAR w/ BART")
print(waic)

## ----warpDLM------------------------------------------------------------------
#Visualize the data
plot(discoveries)

#Fit the model
warpfit <- warpDLM(y=discoveries, type="trend")

## ----warpPPC------------------------------------------------------------------
ppc_ribbon(y=as.vector(discoveries), yrep=warpfit$post_pred)

